
import modefilters
import numpy as np
import pyfits
import read_params
import os,sys,fnmatch,inspect
# import matplotlib
# matplotlib.use(u'Agg')
from matplotlib import pyplot as plt
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject, Gdk
from scipy import signal,linalg,fftpack,integrate,interpolate
import fnmatch
import plotc

Rsun = 695.8
Lx = read_params.get_xlength()
nx = read_params.get_nx()
dt_sec = read_params.get_dt()
dt_min = dt_sec/60
datadir = read_params.get_directory()
codedir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
masterpixels = np.loadtxt(os.path.join(datadir,"master.pixels"),ndmin=1)
solartime = read_params.get_total_simulation_time()

def get_resource_path(rel_path):
    dir_of_py_file = os.path.dirname(inspect.getsourcefile(lambda:0))
    rel_path_to_resource = os.path.join(dir_of_py_file, rel_path)
    abs_path_to_resource = os.path.abspath(rel_path_to_resource)
    return abs_path_to_resource

class MainWindow(Gtk.Window):

    def __init__(self):

        # Files to be used to compute travel times
        self.data_file = None
        self.vzcc_file = None
        self.mode_filter = None
        self.mode_speeds = [0.5,0.75,0.95,1.15,1.2,1.4,1.7,1.9]
        self.source_x = 0
        self.nt = 40+int(solartime*60/dt_min)
        self.t = np.arange(self.nt)*dt_sec
        self.x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)
        self.freq_mHz = np.fft.fftfreq(self.nt,dt_sec)*1e3
        self.k_invMm = 2*np.pi*np.fft.fftfreq(nx,Lx/nx)
        self.halftime = 90 # Minutes

        Gtk.Window.__init__(self, title="Traveltimes")
        self.set_size_request(400, 200)
        self.set_resizable(False)
        try:
            self.set_icon_from_file(get_resource_path("sun.ico"))
        except: pass

        self.timeout_id = None

        vbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=5)
        self.add(vbox)

        align_list = []
        hbox_list = []

        def align_and_hbox():
            align_list.append(Gtk.Alignment(xalign=0.2, yalign=0.0, xscale=0.0, yscale=0.0))
            vbox.add(align_list[-1])
            hbox_list.append(Gtk.Box(spacing=10))
            align_list[-1].add(hbox_list[-1])

        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        self.src_label = Gtk.Label("Source")
        hbox_list[-1].pack_start(self.src_label,expand=True,fill=False,padding=0)

        adjustment_src = Gtk.Adjustment(value=abs(masterpixels).argmin()+1, lower=1, upper=8,
                                        step_incr=1, page_incr=0, page_size=0)
        self.src = Gtk.SpinButton()
        self.src.set_adjustment(adjustment_src)
        self.src.set_numeric(True)
        self.src.set_snap_to_ticks(True)
        hbox_list[-1].pack_start(self.src,expand=True,fill=False,padding=0)
        self.src.connect("value-changed",self.on_src_changed)

        self.source_x = masterpixels[self.src.get_value_as_int()-1]

        self.src_position_label = Gtk.Label("Source located at {:2.0f} Mm".format(self.source_x))
        hbox_list[-1].pack_start(self.src_position_label,expand=True,fill=False,padding=0)

        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        self.iter_label = Gtk.Label("Iter no")
        hbox_list[-1].pack_start(self.iter_label,expand=True,fill=False,padding=0)

        # Get total number of available iterations
        total_iters = len(fnmatch.filter(os.listdir(os.path.join(datadir,"tt")),"iter[0-9][0-9]"))

        adjustment_iter = Gtk.Adjustment(value=0, lower=0, upper=total_iters -1,
                                        step_incr=1, page_incr=0, page_size=0)
        self.iter_no = Gtk.SpinButton()
        self.iter_no.set_adjustment(adjustment_iter)
        self.iter_no.set_numeric(True)
        self.iter_no.set_snap_to_ticks(True)
        self.iter_no.connect("value-changed",self.on_iter_changed)
        hbox_list[-1].pack_start(self.iter_no,expand=True,fill=False,padding=0)

        self.total_iter_no_label = Gtk.Label("Total no of iterations: {:2d}".format(total_iters))
        hbox_list[-1].pack_start(self.total_iter_no_label,expand=True,fill=False,padding=0)

        self.get_data_and_vzcc()



        # Load filter
        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        ridge_filters = ["f"]
        ridge_filters.extend(["p"+str(i) for i in range(1,8)])
        self.ridge_filters_list = Gtk.ComboBoxText()
        self.ridge_filters_list.set_entry_text_column(0)
        for ridge in ridge_filters: self.ridge_filters_list.append_text(ridge)
        self.ridge_filters_list.set_active(0)
        self.ridge_filters_list.connect("changed",self.on_filter_change)

        self.filter_label = Gtk.Label("Filter")
        hbox_list[-1].pack_start(self.filter_label,expand=False,fill=False,padding=0)

        hbox_list[-1].pack_start(self.ridge_filters_list,expand=False,fill=False,padding=10)

        self.params_dist = {"min":20,"max":100}
        self.params_label = Gtk.Label(" "*60)
        # self.on_filter_change(self.ridge_filters_list)
        hbox_list[-1].pack_start(self.params_label,expand=False,fill=False,padding=10)

        # Select pixel

        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        self.pixel_label = Gtk.Label("Pixel:  ")
        hbox_list[-1].pack_start(self.pixel_label,expand=False,fill=False,padding=0)

        starting_x = masterpixels[self.src.get_value_as_int()-1]+20
        starting_pixel = int((starting_x/Lx+0.5)*nx)

        adjustment_pix = Gtk.Adjustment(value=starting_pixel, lower=0, upper=nx-1,
                                        step_incr=1, page_incr=0, page_size=0)
        self.pixel = Gtk.SpinButton()
        self.pixel.set_adjustment(adjustment_pix)
        self.pixel.set_numeric(True)
        self.pixel.set_snap_to_ticks(True)
        self.pixel.connect("value-changed",self.on_pixel_changed)
        hbox_list[-1].pack_start(self.pixel,expand=False,fill=False,padding=0)

        self.x_label = Gtk.Label("x (Mm):")
        hbox_list[-1].pack_start(self.x_label,expand=False,fill=False,padding=0)

        adjustment_x = Gtk.Adjustment(value=starting_x, lower=self.x[0], upper=self.x[-1],
                                        step_incr=self.x[1]-self.x[0], page_incr=0, page_size=0)
        self.xcoord = Gtk.SpinButton()
        self.xcoord.set_adjustment(adjustment_x)
        self.xcoord.set_numeric(True)
        self.xcoord.set_snap_to_ticks(True)
        self.xcoord.connect("value-changed",self.on_xcoord_changed)
        hbox_list[-1].pack_start(self.xcoord,expand=False,fill=False,padding=0)

        # Distance label
        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        self.distance_label = Gtk.Label("Dist from source: {:.1f} Mm,"
                            " min : {:.1f}, max: {:.1f}".format(
                            abs(self.xcoord.get_value() - self.source_x),
                            self.params_dist["min"],self.params_dist["max"] ))
        self.on_filter_change(self.ridge_filters_list)
        hbox_list[-1].pack_start(self.distance_label,expand=False,fill=False,padding=0)


        # Compute button
        align_and_hbox()
        self.pixel_tt_button = Gtk.Button.new_with_label("Pixel Travel Time")
        hbox_list[-1].pack_start(self.pixel_tt_button,expand=True,fill=True,padding=0)
        self.pixel_tt_button.connect("clicked",self.on_pixel_tt_pressed)

        self.all_pixels_tt_button = Gtk.Button.new_with_label("All pixels travel time")
        hbox_list[-1].pack_start(self.all_pixels_tt_button,expand=True,fill=True,padding=0)
        self.all_pixels_tt_button.connect("clicked",self.on_compute_all_pressed)

    def on_pixel_changed(self,pixel):
        # Check if within params range
        current_pixel = pixel.get_value()
        current_x = (current_pixel/nx-0.5)*Lx
        self.xcoord.set_value(current_x)
        self.distance_label.set_text("Dist from source: {:.1f} Mm,"
                            " min : {:.1f}, max: {:.1f}".format(
                            abs(self.xcoord.get_value() - self.source_x),
                            self.params_dist["min"],self.params_dist["max"]))

    def on_xcoord_changed(self,xcoord):
        # Check if within params range

        current_x = xcoord.get_value()
        current_pixel = int(nx//2 + current_x/Lx*nx)
        self.pixel.set_value(current_pixel)
        self.distance_label.set_text("Dist from source: {:.1f} Mm,"
                            " min : {:.1f}, max: {:.1f}".format(
                            abs(self.xcoord.get_value() - self.source_x),
                            self.params_dist["min"],self.params_dist["max"]))

    def on_filter_change(self,ridge_filters_list):
        # Try to load params file to get min and max dist
        active_index = ridge_filters_list.get_active()
        params_path = os.path.join(datadir,"params."+str(active_index))
        try:
            params = np.loadtxt(params_path,ndmin=1)
            self.params_dist["min"] = params[0]
            self.params_dist["max"] = params[1]
            self.halftime = int(params[2]/(2*dt_min))

        except IOError:
            self.params_dist["min"] = 20
            self.params_dist["max"] = 100
            self.halftime = 90

        # self.params_label.set_text("Min dist: {:.1f}, max dist: {:.1f}".
        #                             format(self.params_dist["min"],self.params_dist["max"]))
        self.distance_label.set_text("Dist from source: {:.1f} Mm,"
                            " min : {:.1f}, max: {:.1f}".format(
                            abs(self.xcoord.get_value() - self.source_x),
                            self.params_dist["min"],self.params_dist["max"]))
        filter_file_path = os.path.join(codedir,
                    self.ridge_filters_list.get_active_text()+"mode_filter.fits")
        try:
            self.mode_filter = np.squeeze(pyfits.getdata(filter_file_path)).astype(float)
        except IOError:
            modefilter_function = getattr(modefilters,
                    self.ridge_filters_list.get_active_text()+"mode_filter")
            self.mode_filter = np.squeeze(modefilter_function(self.nt,dt_sec,nx,Lx)).T

    def on_data_entry_changed(self,entry):
        current_path = entry.get_text()
        if os.path.isfile(current_path) and os.path.splitext(current_path)[1].lower()==".fits":
            try:
                self.data_file = np.squeeze(pyfits.getdata(current_path)).astype(float)
            except IOError:
                self.data_file = None

    def on_vzcc_entry_changed(self,entry):
        current_path = entry.get_text()
        if os.path.isfile(current_path) and os.path.splitext(current_path)[1].lower()==".fits":
            try:
                self.vzcc_file = np.squeeze(pyfits.getdata(current_path)).astype(float)
            except IOError:
                self.vzcc_file = None

    def on_src_changed(self,src):
        self.source_x = masterpixels[self.src.get_value_as_int()-1]
        self.get_data_and_vzcc()
        self.src_position_label.set_text("Source located at {:2.0f} Mm".format(self.source_x))

    def on_iter_changed(self,widget):
        self.get_data_and_vzcc(data=False)

    def get_data_and_vzcc(self,data=True):
        # Load data
        if data:
            try:
                src_file_path = os.path.join(datadir,"data",
                                str(self.src.get_value_as_int()).zfill(2)+".fits")
                self.data_file = fitsread(src_file_path)
                self.nt = self.data_file.shape[0]
                self.time_coordinates = np.arange(self.nt)*dt_sec
            except IOError:
                self.data_file = None

        # Load vzcc
        try:
            vzcc_file_path = os.path.join(datadir,"tt",
                            "iter"+str(self.iter_no.get_value_as_int()).zfill(2),
                            "vz_cc_src"+str(self.src.get_value_as_int()).zfill(2)+".fits")
            self.vzcc_file = fitsread(vzcc_file_path)

        except IOError:
            self.vzcc_file = None

    def filter_data_vzcc(self):
        if self.mode_filter is not None:
            data_filtered = np.fft.ifft2(np.fft.fft2(self.data_file)*self.mode_filter).real
            vzcc_filtered = np.fft.ifft2(np.fft.fft2(self.vzcc_file)*self.mode_filter).real
        else:
            data_filtered = self.data_file
            vzcc_filtered = self.vzcc_file
        return data_filtered,vzcc_filtered

    def on_pixel_tt_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None:
            return

        pixel = self.pixel.get_value_as_int()

        data_filtered,vzcc_filtered = self.filter_data_vzcc()

        distance_Mm = abs(self.xcoord.get_value() - self.source_x)
        signed_distance_Mm = self.xcoord.get_value() - self.source_x
        if distance_Mm<self.params_dist["min"]:
            print("Distance",distance_Mm,"is too low to compute travel time")
            return
        elif distance_Mm>self.params_dist["max"]:
            print("Distance",distance_Mm,"is too high to compute travel time")
            return

        print(pixel,"Distance from source",signed_distance_Mm,"Mm")
        u0 = data_filtered[:,pixel]
        u = vzcc_filtered[:,pixel]

        vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
        timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))
        timefin_index = timest_index + 60

        loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index
        t_low_index = max(0,loc - self.halftime)
        t_high_index = min(self.nt,loc + 2*self.halftime)

        tt_han = self.compute_tt_quadratic_fit(u0,u,t_low_index,t_high_index)
        print("Hanasoge",tt_han,"seconds")
        tt_gb = self.compute_tt_gizonbirch(u0,u,t_low_index,t_high_index)
        print("Gizon Birch",tt_gb,"seconds")

        plt.figure(0)
        plt.plot(self.t/60,u0,label="data",color="blue")
        plt.plot(self.t/60,u,label="vzcc",color="red")
        # plt.axvline(timest_index*dt_min,ls="dotted",color="black")
        # plt.axvline(timefin_index*dt_min,ls="dotted",color="black")
        plt.axvspan(t_low_index*dt_min,t_high_index*dt_min,edgecolor="0.6",facecolor="0.95")
        plt.xlim(self.t[0]/60,self.t[-1]/60)
        plt.xlabel("Time (min)",fontsize=16)
        plt.ylabel("Wave displacement",fontsize=16)
        plt.legend(loc="best",fontsize=14)
        plt.tight_layout()

        plt.figure(1)
        plotc.colorplot(trig_interp(data_filtered,x=4,y=2),
        x=lininterp(self.x,factor=4),y=lininterp(self.t,factor=2)/60,
        sp=121,colorbar=False)
        plt.axvline(self.x[pixel],linestyle="dashed",linewidth=1,color="firebrick")
        plt.axhline(t_low_index*dt_min,color="0.6")
        plt.axhline(t_high_index*dt_min,color="0.6")
        plt.xlim(max(-Lx/2,self.x[pixel]-15),min(Lx/2,self.x[pixel]+15))
        plt.ylim(t_low_index*dt_min-10,t_high_index*dt_min+10)
        plt.title("Time-Distance",fontsize=16)
        plt.xlabel("x (Mm)",fontsize=14)
        plt.ylabel("Time (min)",fontsize=14)

        plotc.spectrumplot(abs(np.fft.fft2(data_filtered)),x=self.k_invMm*Rsun,y=self.freq_mHz,
                    sp=122,cmap="Greys",xr=[0,1200],yr=[0,8],colorbar=False)
        plt.title("Spectrum",fontsize=16)
        plt.xlabel(r"$k\,R_\odot$",fontsize=14)
        plt.ylabel(r"$\nu$ (mHz)",fontsize=14)

        plt.tight_layout()
        plt.show()

    def on_compute_all_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None: return

        data_filtered,vzcc_filtered = self.filter_data_vzcc()

        tt_hanasoge_30s = []
        tt_gb_30s = []
        pixel_x_list = []

        for pixel in range(nx):

            pixel_x = (pixel/nx-0.5)*Lx
            distance_Mm = abs(pixel_x - self.source_x)
            signed_distance_Mm = pixel_x - self.source_x
            if distance_Mm<self.params_dist["min"] or distance_Mm>self.params_dist["max"]:
                continue

            vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
            timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))-1
            timefin_index = timest_index + 60

            if timest_index>self.nt:
                print("Pixel {:3d}, starting index {:3d} is greater than nt {:3d}"
                .format(pixel,timest_index,self.nt))
                continue
            if timefin_index>self.nt:
                print("Pixel {:3d}, ending index {:3d} is greater than nt {:3d}\nChoosing nt instead"
                        .format(pixel,timest_index,self.nt))
                timefin_index = self.nt

            u0 = data_filtered[:,pixel]
            u = vzcc_filtered[:,pixel]

            loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index

            t_low_index = max(0,loc - self.halftime)
            t_high_index = min(self.nt,loc + 2*self.halftime)

            tt = self.compute_tt_quadratic_fit(u0,u,t_low_index,t_high_index)
            tt_hanasoge_30s.append(tt)

            tt = self.compute_tt_gizonbirch(u0,u,t_low_index,t_high_index)
            tt_gb_30s.append(tt)
            pixel_x_list.append(signed_distance_Mm)

        pixel_x_list = np.array(pixel_x_list)
        tt_gb_30s = np.array(tt_gb_30s)
        tt_hanasoge_30s = np.array(tt_hanasoge_30s)

        color = iter(['steelblue','olive'])

        plt.cla()

        def plot_tt(arr,**kwargs):
            left_region = (pixel_x_list>-self.params_dist["max"]) & (pixel_x_list<-self.params_dist["min"])
            p = plt.plot(pixel_x_list[left_region],arr[left_region],color=kwargs.pop("color",next(color)),
            marker=kwargs.pop("marker",'.'),**kwargs)

            right_region = (pixel_x_list>self.params_dist["min"]) & (pixel_x_list<self.params_dist["max"])
            plt.plot(pixel_x_list[right_region],arr[right_region],color=p[0].get_color(),
            marker=kwargs.pop("marker",p[0].get_marker()))

        plot_tt(tt_hanasoge_30s,label="Sparc, 30s")
        plot_tt(tt_gb_30s,label="GB02, 30s")

        plt.axvline(self.source_x,color="firebrick",label="Source")
        plt.axvspan(-15,15,facecolor="0.95",edgecolor="0.6",label="SG")

        plt.margins(x=0.1,y=0.1)
        plt.legend(loc="best",fontsize=14)
        plt.xlabel("Distance from source (Mm)",fontsize=16)
        plt.ylabel(r"$\Delta \tau$ (sec)",fontsize=16)
        plt.title("Source at x = {:2.1f} Mm".format(self.source_x))
        plt.show()

    def compute_tt_quadratic_fit(self,u0, u,t_low_index,t_high_index):
        u0 = np.take(u0,list(range(t_low_index,t_high_index+1)),mode='wrap')
        u  =  np.take(u,list(range(t_low_index,t_high_index+1)),mode='wrap')
        cc=signal.correlate(u0,u).real
        cc_max_index = cc.argmax()
        t = np.arange(-len(u)+1,len(u))*dt_sec

        t_wavepacket = t[cc_max_index-1:cc_max_index+2]
        cc_wavepacket = cc[cc_max_index-1:cc_max_index+2]

        p = np.polyfit(t_wavepacket,cc_wavepacket,2)

        return p[1]/(2*p[0])

    def compute_tt_gizonbirch(self,u0,u,t_low_index,t_high_index):

        w = np.zeros(self.nt,dtype=float)
        np.put(w,list(range(t_low_index,t_high_index+1)),[1],mode='wrap')
        # w[t_low_index:t_high_index+1] = 1

        degree = 4

        u0dot = np.array([fftpack.diff(u0,period=self.nt*dt_sec,order=i) for i in range(degree+2)])

        p = np.zeros(degree)

        p[-1] = integrate.simps(w*2*u0dot[1]*(u-u0),dx=dt_sec)
        p[-2] = integrate.simps(w*2*(u0dot[1]**2+u0dot[2]*(u0-u)),dx=dt_sec)
        p[-3] = integrate.simps(w*(-3*u0dot[1]*u0dot[2]+u0dot[3]*(u-u0)),dx=dt_sec)
        p[-4] = integrate.simps(
                w*(3*u0dot[2]**2+4*u0dot[1]*u0dot[3]+u0dot[4]*(u0-u))/3,dx=dt_sec)

        guess = -p[-1]/p[-2]
        r = np.roots(p)

        return r[abs(r-guess).argmin()]


def close_plot_and_exit(*args):
    plt.close('all')
    Gtk.main_quit(args)

def fitsread(f): return np.squeeze(pyfits.getdata(f)).astype(float)

def trig_interp(arr,x=1,y=1):
    if arr.ndim==1:
        spec = np.fft.fft(arr)
        nx, = spec.shape
        spec = np.concatenate((spec[:nx//2],
                            np.zeros(nx*(x-1),dtype=complex),
                            spec[nx//2:]))

        return np.fft.ifft(spec).real

    elif arr.ndim==2:
        spec = np.fft.fft2(arr)
        ny,nx = spec.shape

        top = np.concatenate((spec[:ny//2,:nx//2],
                            np.zeros((ny//2,nx*(x-1)),dtype=complex),
                            spec[:ny//2,nx//2:]),axis=1)

        middle = np.zeros((ny*(y-1),nx*x),dtype=complex)

        bottom = np.concatenate((spec[ny//2:,:nx//2],
                            np.zeros((ny//2,nx*(x-1)),dtype=complex),
                            spec[ny//2:,nx//2:]),axis=1)

        spec = np.concatenate((top,middle,bottom),axis=0)

        return np.fft.ifft2(spec).real

def lininterp(arr,factor=1):
    dx = arr[1]-arr[0]
    nx, = arr.shape
    return np.linspace(arr[0],arr[-1]+dx,nx*factor,endpoint=False)

win1 = MainWindow()
win1.connect("delete-event", close_plot_and_exit)
win1.connect('destroy', Gtk.main_quit)
win1.show_all()
Gtk.main()
