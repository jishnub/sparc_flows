from __future__ import division,print_function
import modefilters
import numpy as np
import pyfits
import read_params
import os,sys,fnmatch,inspect
from matplotlib import pyplot as plt
import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject, Gdk
from scipy import signal,linalg,fftpack,integrate,interpolate
import fnmatch

Rsun = 695.8
Lx = read_params.get_xlength()
nx = read_params.get_nx()
dt_sec = read_params.get_dt()
dt_min = dt_sec/60
nu_max = 1e3/dt_sec/2
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
        self.freq_filter = None
        self.mode_speeds = [0.5,0.75,0.95,1.15,1.2,1.4,1.7,1.9] # Mm/min
        self.mode_freq_ranges = {"f":(2,3),
        "p1":(2.7,4.2),"p2":(3.1,5.2),"p3":(3.4,5.8),"p4":(3,6),"p5":(3.5,6),"p6":(3.5,6),"p7":(3.5,6)}
        self.source_x = 0
        self.nt = 40+int(solartime*60/dt_min)
        self.t = np.arange(self.nt)*dt_sec
        self.halftime = 10 # Minutes
        self.x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

        Gtk.Window.__init__(self, title="Traveltimes")
        self.set_size_request(400, 230)
        self.set_resizable(False)
        self.set_icon_from_file(get_resource_path("sun.ico"))

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
        ridge_filters.extend(["p"+str(i) for i in xrange(1,8)])
        self.ridge_filters_list = Gtk.ComboBoxText()
        self.ridge_filters_list.set_entry_text_column(0)
        for ridge in ridge_filters: self.ridge_filters_list.append_text(ridge)

        # Get a list of params files in the data directory
        params=fnmatch.filter(os.listdir(read_params.get_directory()),"params.[0-9]")
        self.ridge_filters_list.set_active(int(os.path.splitext(params[0])[1][1:]))
        self.ridge_filters_list.connect("changed",self.on_filter_change)

        self.filter_label = Gtk.Label("Filter")
        hbox_list[-1].pack_start(self.filter_label,expand=False,fill=False,padding=0)
        hbox_list[-1].pack_start(self.ridge_filters_list,expand=False,fill=False,padding=0)

        # Freq filter
        self.freq_filter_label = Gtk.Label("Freq")
        hbox_list[-1].pack_start(self.freq_filter_label,expand=False,fill=False,padding=0)

        freq_low_iter = Gtk.Adjustment(
        value=self.mode_freq_ranges[self.ridge_filters_list.get_active_text()][0],
        lower=2.0, upper=6.5,step_incr=0.5)
        self.freq_low = Gtk.SpinButton()
        self.freq_low.set_adjustment(freq_low_iter)
        self.freq_low.set_digits(1)
        self.freq_low.set_numeric(True)
        self.freq_low.set_snap_to_ticks(True)
        self.freq_low.connect("value-changed",self.on_freq_low_changed)
        hbox_list[-1].pack_start(self.freq_low,expand=False,fill=False,padding=0)

        freq_high_iter = Gtk.Adjustment(
        value=self.mode_freq_ranges[self.ridge_filters_list.get_active_text()][1],
        lower=2.0, upper=6.5,step_incr=0.5)
        self.freq_high = Gtk.SpinButton()
        self.freq_high.set_adjustment(freq_high_iter)
        self.freq_high.set_digits(1)
        self.freq_high.set_numeric(True)
        self.freq_high.set_snap_to_ticks(True)
        self.freq_high.connect("value-changed",self.on_freq_high_changed)
        hbox_list[-1].pack_start(self.freq_high,expand=False,fill=False,padding=0)

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

        self.tt_label = Gtk.Label("TT window")
        hbox_list[-1].pack_start(self.tt_label, False, False, 0)

        self.tt_velocity = Gtk.RadioButton.new_with_label_from_widget(None, "Mode Velocity")
        hbox_list[-1].pack_start(self.tt_velocity, False, False, 0)

        self.tt_peak_locate = Gtk.RadioButton.new_from_widget(self.tt_velocity)
        self.tt_peak_locate.set_label("Locate Peak")
        hbox_list[-1].pack_start(self.tt_peak_locate, False, False, 0)



        # Distance label
        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        self.params_dist = {"min":20,"max":100}
        self.distance_label = Gtk.Label()
        self.set_distance_label()
        self.on_filter_change(self.ridge_filters_list)
        hbox_list[-1].pack_start(self.distance_label,expand=False,fill=False,padding=0)


        # Compute button
        align_and_hbox()
        self.pixel_tt_button = Gtk.Button.new_with_label("Single Pixel")
        hbox_list[-1].pack_start(self.pixel_tt_button,expand=True,fill=True,padding=0)
        self.pixel_tt_button.connect("clicked",self.on_pixel_tt_pressed)

        self.all_pixels_tt_button = Gtk.Button.new_with_label("All pixels")
        hbox_list[-1].pack_start(self.all_pixels_tt_button,expand=True,fill=True,padding=0)
        self.all_pixels_tt_button.connect("clicked",self.on_compute_all_pressed)

        self.freq_montage_button = Gtk.Button.new_with_label("Freq Montage")
        hbox_list[-1].pack_start(self.freq_montage_button,expand=True,fill=True,padding=0)
        self.freq_montage_button.connect("clicked",self.on_freq_montage_pressed)

    def on_pixel_changed(self,pixel):
        # Check if within params range
        current_pixel = pixel.get_value()
        current_x = (current_pixel/nx-0.5)*Lx
        self.xcoord.set_value(current_x)
        self.set_distance_label()

    def on_xcoord_changed(self,xcoord):
        # Check if within params range

        current_x = xcoord.get_value()
        current_pixel = int(nx//2 + current_x/Lx*nx)
        self.pixel.set_value(current_pixel)
        self.set_distance_label()

    def on_filter_change(self,ridge_filters_list):
        # Try to load params file to get min and max dist
        active_index = ridge_filters_list.get_active()
        filter_name = ridge_filters_list.get_active_text()
        params_path = os.path.join(datadir,"params."+str(active_index))
        try:
            params = np.loadtxt(params_path,ndmin=1)
            self.params_dist["min"] = params[0]
            self.params_dist["max"] = params[1]
            self.halftime = int(params[2]/(2*dt_min))

        except IOError:
            self.params_dist["min"] = 20
            self.params_dist["max"] = 100
            self.halftime = 10

        self.set_distance_label()

        modefilter_function = getattr(modefilters,
                self.ridge_filters_list.get_active_text()+"mode_filter")
        self.mode_filter = np.squeeze(modefilter_function(self.nt,dt_sec,nx,Lx)).T

        self.freq_low.set_value(self.mode_freq_ranges[filter_name][0])
        self.freq_high.set_value(self.mode_freq_ranges[filter_name][1])


    def on_freq_low_changed(self,SpinButton):
        step_size,_ = self.freq_high.get_increments()
        if self.freq_high.get_value()<self.freq_low.get_value()+step_size:
            self.freq_high.set_value(self.freq_low.get_value()+step_size)

    def on_freq_high_changed(self,SpinButton):
        step_size,_ = self.freq_low.get_increments()
        if self.freq_low.get_value()>self.freq_high.get_value()-step_size:
            self.freq_low.set_value(self.freq_high.get_value()-step_size)

    def get_freq_filter(self,low=None,high=None):
        if low is None: low=self.freq_low.get_value()
        if high is None: high=self.freq_high.get_value()
        return modefilters.freq_filter(low,high,self.nt,dt_sec)[:,None]

    def get_time_window(self,u):
        max_loc = abs(u).argmax()
        # find peaks around maximum

        window = 40
        low_t = max(max_loc-window,0)
        high_t = max_loc + window
        u_window = u[low_t:high_t]
        peak_locations = find_maxima(u_window)+1 #+low_t+1 # Indices of u instead of u_window
        peaks = np.log(u_window[peak_locations])
        p = np.polyfit(peak_locations+low_t,peaks,2)

        x0 = -p[1]/(2*p[0])
        sigma = 1/np.sqrt(-2*p[0])
        # A = np.exp(p[2]+x0**2/(2*sigma**2))
        arrival = max(int(x0-sigma),0)
        departure = min(int(x0+sigma),self.nt-1)
        return arrival,departure


    def set_distance_label(self):
        self.distance_label.set_text("Dist from source: {:.1f} Mm,"
                            " min : {:.1f}, max: {:.1f}".format(
                            abs(self.xcoord.get_value() - self.source_x),
                            self.params_dist["min"],self.params_dist["max"]))

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
                self.t = np.arange(self.nt)*dt_sec
                self.freq_filter = np.ones(self.nt)
                print("Loaded data file",src_file_path)
            except IOError:
                print("Could not load data file")
                self.data_file = None

        # Load vzcc
        try:
            vzcc_file_path = os.path.join(datadir,"tt",
                            "iter"+str(self.iter_no.get_value_as_int()).zfill(2),
                            "vz_cc_src"+str(self.src.get_value_as_int()).zfill(2)+".fits")
            self.vzcc_file = fitsread(vzcc_file_path)
            print("Loaded vzcc file",vzcc_file_path)

        except IOError:
            print("Could not load vzcc file")
            self.vzcc_file = None

        try:
            self.set_distance_label()
        except AttributeError: pass

    def filter_data_vzcc(self,low=None,high=None):
        if low is None: low=self.freq_low.get_value()
        if high is None: high=self.freq_high.get_value()
        self.freq_filter = self.get_freq_filter(low,high)
        if self.mode_filter is not None:
            data_filtered = np.fft.ifft2(np.fft.fft2(self.data_file)*self.mode_filter*self.freq_filter).real
            vzcc_filtered = np.fft.ifft2(np.fft.fft2(self.vzcc_file)*self.mode_filter*self.freq_filter).real
        else:
            data_filtered = self.data_file
            vzcc_filtered = self.vzcc_file
        return data_filtered,vzcc_filtered

    def on_pixel_tt_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None:
            print("No data or vzcc loaded")
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

        print()
        print("Src {:1d}, iter {:2d}, pixel {:3d}, distance from source {:.1f} Mm".format(
        self.src.get_value_as_int(),self.iter_no.get_value_as_int(),pixel,signed_distance_Mm))
        print(self.ridge_filters_list.get_active_text(),"mode, freq filter",
        self.freq_low.get_value(),"to",self.freq_high.get_value())
        u0 = data_filtered[:,pixel]
        u = vzcc_filtered[:,pixel]

        if self.tt_velocity.get_active():
            vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
            timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))
            timefin_index = timest_index + 40

            loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index
            t_low_index = loc - self.halftime
            t_high_index = loc + self.halftime
        elif self.tt_peak_locate.get_active():
            timest_index,timefin_index,loc = None,None,None
            t_low_index,t_high_index = self.get_time_window(u0)

        tt_han = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=1)
        print("Hanasoge",tt_han,"seconds")
        tt_han_smoothed = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=3)
        print("Smoothed, 10s cadence",tt_han_smoothed,"seconds")
        # tt_han_smoothed = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=6)
        # print("Smoothed, 5s cadence",tt_han_smoothed,"seconds")
        tt_gb = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=1)
        print("Gizon Birch",tt_gb,"seconds")
        # tt_gb_interp = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=3)
        # print("Gizon Birch, 10s",tt_gb_interp,"seconds")
        # tt_gb_interp = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=6)
        # print("Gizon Birch, 5s",tt_gb_interp,"seconds")

        print("Fortran array indices to compare with sparc output")
        print(("{:6s}{:6s}"+"{:6s}"*5).format("pixel","tt","lef","rig","loc","tst","tfin"))
        print(("{:<6d}{:<6.1f}"+"{:<6d}"*5).format(
        pixel+1,tt_gb,t_low_index+1,t_high_index+1,loc+1,timest_index,timefin_index))

        def handle_close(evt):
            plt.close('all')

        tdfig = plt.figure()
        tdfig.canvas.mpl_connect('close_event', handle_close)
        tdfig_ax = tdfig.gca()
        p=tdfig_ax.pcolormesh(self.x,self.t/3600,data_filtered,cmap="RdBu_r")
        tdfig_ax.axhline(t_low_index*dt_min/60,ls="dotted",color="black")
        tdfig_ax.axhline(t_high_index*dt_min/60,ls="dotted",color="black")
        tdfig_ax.axvline(self.xcoord.get_value(),ls="dashed",color="black")
        tdfig_ax.set_xlabel("x (Mm)")
        tdfig_ax.set_ylabel("time (Hours)")
        tdfig_ax.set_xlim(self.x[0],self.x[-1])
        tdfig_ax.set_ylim(self.t[0]/3600,self.t[-1]/3600)
        tdfig_ax.set_title(self.ridge_filters_list.get_active_text()+" mode time distance")
        tdfig.colorbar(p,format="%.1e")

        ttfig = plt.figure()
        ttfig.canvas.mpl_connect('close_event', handle_close)
        ttfig_ax = ttfig.gca()

        ttfig_ax.plot(self.t/60,u0,label="data",color="blue")
        ttfig_ax.plot(self.t/60,u,label="vzcc",color="red")
        if self.tt_velocity.get_active():
            ttfig_ax.axvline(timest_index*dt_min,ls="dashed",color="black")
            ttfig_ax.axvline(timefin_index*dt_min,ls="dashed",color="black")
        ttfig_ax.axvspan(t_low_index*dt_min,t_high_index*dt_min,color="0.8")
        ttfig_ax.set_xlim(0,self.nt*dt_min)
        ttfig_ax.set_title("Pixel: {:3d}, dist from src: {:.1f} Mm, {:s} mode"
        r" $\Delta\tau$ (GB02): {:.1f} sec".format(pixel,distance_Mm,
        self.ridge_filters_list.get_active_text(),tt_gb),fontsize=14)
        ttfig_ax.set_xlabel("Time (min)",fontsize=16)
        ttfig_ax.set_ylabel("Wave displacement",fontsize=16)
        ttfig_ax.legend(loc="best",fontsize=14)

        plt.tight_layout()
        plt.show()

    def on_compute_all_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None:
            print("No data or vzcc loaded")
            return

        data_filtered,vzcc_filtered = self.filter_data_vzcc()

        tt_hanasoge_30s = []
        tt_gb_30s = []
        tt_gb_interp_10s = []
        tt_gb_interp_5s = []
        tt_hanasoge_interp_10s = []
        tt_hanasoge_interp_5s = []
        pixel_x_list = []

        for pixel in xrange(nx):

            pixel_x = (pixel/nx-0.5)*Lx
            distance_Mm = abs(pixel_x - self.source_x)
            signed_distance_Mm = pixel_x - self.source_x
            if distance_Mm<self.params_dist["min"] or distance_Mm>self.params_dist["max"]:
                continue

            u0 = data_filtered[:,pixel]
            u = vzcc_filtered[:,pixel]

            vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
            timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))-1
            timefin_index = timest_index + 40

            if self.tt_velocity.get_active():
                if timest_index>self.nt:
                    # print("Pixel {:3d}, starting index {:3d} is greater than nt {:3d}"
                    # .format(pixel,timest_index,self.nt))
                    continue
                if timefin_index>self.nt:
                    # print("Pixel {:3d}, ending index {:3d} is greater than nt {:3d}\nChoosing nt instead"
                    #         .format(pixel,timest_index,self.nt-1))
                    timefin_index = self.nt-1
                loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index

                t_low_index = loc - self.halftime
                t_high_index = loc + self.halftime

            elif self.tt_peak_locate.get_active():
                t_low_index,t_high_index = self.get_time_window(u0)

            t_low_index = max(0,t_low_index)
            t_high_index = min(self.nt-1,t_high_index)

            tt = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=1)
            tt_hanasoge_30s.append(tt)
            tt = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=3)
            tt_hanasoge_interp_10s.append(tt)
            tt = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=6)
            tt_hanasoge_interp_5s.append(tt)
            tt = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=1)
            tt_gb_30s.append(tt)
            tt = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=3)
            tt_gb_interp_10s.append(tt)
            tt = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=6)
            tt_gb_interp_5s.append(tt)
            pixel_x_list.append(signed_distance_Mm)

        pixel_x_list = np.array(pixel_x_list)
        tt_gb_30s = np.array(tt_gb_30s)
        tt_gb_interp_10s = np.array(tt_gb_interp_10s)
        tt_gb_interp_5s = np.array(tt_gb_interp_5s)
        tt_hanasoge_30s = np.array(tt_hanasoge_30s)
        tt_hanasoge_interp_10s = np.array(tt_hanasoge_interp_10s)
        tt_hanasoge_interp_5s = np.array(tt_hanasoge_interp_5s)

        # tt_file=np.loadtxt(os.path.join(datadir,"forward_src"+str(self.src.get_value_as_int()).zfill(2)+
        #                     "_ls00","ttdiff."+str(self.ridge_filters_list.get_active())))

        color = iter(['r','maroon','indianred','b','steelblue','royalblue'])

        plt.cla()

        def plot_tt(arr,**kwargs):
            left_region = (pixel_x_list>-self.params_dist["max"]) & (pixel_x_list<-self.params_dist["min"])
            p = plt.plot(pixel_x_list[left_region],arr[left_region],color=kwargs.pop("color",next(color)),
            marker=kwargs.pop("marker",'.'),**kwargs)

            right_region = (pixel_x_list>self.params_dist["min"]) & (pixel_x_list<self.params_dist["max"])
            plt.plot(pixel_x_list[right_region],arr[right_region],color=p[0].get_color(),
            marker=kwargs.pop("marker",p[0].get_marker()))

        plot_tt(tt_hanasoge_30s,label="Sparc, 30s")
        plot_tt(tt_hanasoge_interp_10s,label="Sparc,interp: 10s")
        plot_tt(tt_hanasoge_interp_5s,label="Sparc,interp: 5s")
        # plot_tt(tt_file[:,1],label="From file",marker="o",ls="None")


        plot_tt(tt_gb_30s,label="GB02, 30s")
        # plot_tt(tt_gb_interp_10s,label="GB02,interp: 10s")
        # plot_tt(tt_gb_interp_5s,label="GB02,interp: 5s")

        plt.legend(loc="best",fontsize=14)
        plt.xlabel("Distance from source (Mm)",fontsize=16)
        plt.ylabel(r"$\Delta \tau$ (sec)",fontsize=16)
        plt.title("Source at x = {:2.1f} Mm, {:s} mode {:.1f} "
        "mHz to {:.1f} mHz".format(self.source_x,self.ridge_filters_list.get_active_text(),
        self.freq_low.get_value(),self.freq_high.get_value()))
        plt.show()

    def on_freq_montage_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None:
            print("No data or vzcc loaded")
            return

        freqs = np.linspace(self.freq_low.get_value(),
        self.freq_high.get_value(),5)
        # freqs = np.linspace(self.mode_freq_ranges[self.ridge_filters_list.get_active_text()][0],
        # self.mode_freq_ranges[self.ridge_filters_list.get_active_text()][1],5)

        spec_fig = plt.figure()
        spec = abs(np.fft.fft2(np.squeeze(self.data_file)))**2
        kR = np.fft.fftfreq(nx,Lx/nx)*2*np.pi*Rsun
        nu = np.fft.fftfreq(self.nt,dt_sec)*1e3
        plt.pcolormesh(np.fft.fftshift(kR),np.fft.fftshift(nu),np.fft.fftshift(spec),
        cmap="Greys",vmax=abs(spec).max()/5)
        plt.xlim(0,1000)
        plt.ylim(1,7)
        plt.xlabel("$kR_\odot$",fontsize=16)
        plt.ylabel("$nu$ (mHz)",fontsize=16)
        for ind,(freqlow,freqhigh) in enumerate(zip(freqs[:-1],freqs[1:])):
            plt.axhspan(freqlow,freqhigh,facecolor="papayawhip",edgecolor="coral",alpha=0.3)

        tt_fig = plt.figure()
        for ind,(freqlow,freqhigh) in enumerate(zip(freqs[:-1],freqs[1:])):
            plt.subplot(2,2,ind+1)
            data_filtered,vzcc_filtered = self.filter_data_vzcc(freqlow,freqhigh)

            tt_hanasoge_30s = []
            tt_gb_30s = []
            pixel_x_list = []

            for pixel in xrange(nx):

                pixel_x = (pixel/nx-0.5)*Lx
                distance_Mm = abs(pixel_x - self.source_x)
                signed_distance_Mm = pixel_x - self.source_x
                if distance_Mm<self.params_dist["min"] or distance_Mm>self.params_dist["max"]:
                    continue

                u0 = data_filtered[:,pixel]
                u = vzcc_filtered[:,pixel]

                vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
                timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))-1
                timefin_index = timest_index + 40

                if self.tt_velocity.get_active():
                    if timest_index>self.nt:
                        # print("Pixel {:3d}, starting index {:3d} is greater than nt {:3d}"
                        # .format(pixel,timest_index,self.nt))
                        continue
                    if timefin_index>self.nt:
                        # print("Pixel {:3d}, ending index {:3d} is greater than nt {:3d}\nChoosing nt instead"
                        #         .format(pixel,timest_index,self.nt-1))
                        timefin_index = self.nt-1
                    loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index

                    t_low_index = loc - self.halftime
                    t_high_index = loc + self.halftime

                elif self.tt_peak_locate.get_active():
                    t_low_index,t_high_index = self.get_time_window(u0)

                t_low_index = max(0,t_low_index)
                t_high_index = min(self.nt-1,t_high_index)

                tt = self.compute_tt_quadratic_fit_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=1)
                tt_hanasoge_30s.append(tt)
                tt = self.compute_tt_gizonbirch_interpolated(u0,u,t_low_index,t_high_index,interpolation_factor=1)
                tt_gb_30s.append(tt)
                pixel_x_list.append(signed_distance_Mm)

            pixel_x_list = np.array(pixel_x_list)
            tt_gb_30s = np.array(tt_gb_30s)
            tt_hanasoge_30s = np.array(tt_hanasoge_30s)

            # tt_file=np.loadtxt(os.path.join(datadir,"forward_src"+str(self.src.get_value_as_int()).zfill(2)+
            #                     "_ls00","ttdiff."+str(self.ridge_filters_list.get_active())))

            color = iter(['peru','teal'])

            plt.cla()

            def plot_tt(arr,**kwargs):
                left_region = (pixel_x_list>-self.params_dist["max"]) & (pixel_x_list<-self.params_dist["min"])
                p = plt.plot(pixel_x_list[left_region],arr[left_region],color=kwargs.pop("color",next(color)),
                marker=kwargs.pop("marker",'.'),**kwargs)

                right_region = (pixel_x_list>self.params_dist["min"]) & (pixel_x_list<self.params_dist["max"])
                plt.plot(pixel_x_list[right_region],arr[right_region],color=p[0].get_color(),
                marker=kwargs.pop("marker",p[0].get_marker()))

            plot_tt(tt_hanasoge_30s,label="Quad fit, 30s")
            plot_tt(tt_gb_30s,label="GB02, 30s")

            plt.legend(loc="best",fontsize=14)
            plt.xlabel("Distance from source (Mm)",fontsize=16)
            plt.ylabel(r"$\Delta \tau$ (sec)",fontsize=16)
            plt.title("{:.1f} mHz to {:.1f} mHz".format(freqlow,freqhigh),fontsize=16)

        plt.tight_layout()

        def handle_close(evt):
            plt.close('all')

        tt_fig.canvas.mpl_connect('close_event', handle_close)
        spec_fig.canvas.mpl_connect('close_event', handle_close)

        plt.show()

    def compute_tt_quadratic_fit_interpolated(self,u0, u,t_low_index,t_high_index,interpolation_factor=1):
        u0 = u0[t_low_index:t_high_index+1]
        u  =  u[t_low_index:t_high_index+1]

        t_fine = np.linspace(self.t[t_low_index],self.t[t_high_index],
                    (t_high_index-t_low_index+1)*interpolation_factor)
        s0 = interpolate.InterpolatedUnivariateSpline(self.t[t_low_index:t_high_index+1],u0)
        u0 = s0(t_fine)
        s = interpolate.InterpolatedUnivariateSpline(self.t[t_low_index:t_high_index+1],u)
        u = s(t_fine)

        cc=signal.correlate(u0,u).real
        cc_max_index = cc.argmax()
        t = np.arange(-len(u)+1,len(u))*dt_sec/interpolation_factor

        t_wavepacket = t[cc_max_index-interpolation_factor:cc_max_index+interpolation_factor+1]
        cc_wavepacket = cc[cc_max_index-interpolation_factor:cc_max_index+interpolation_factor+1]

        p = np.polyfit(t_wavepacket,cc_wavepacket,2)

        return -p[1]/(2*p[0])

    def compute_tt_gizonbirch_interpolated(self,u0,u,t_low_index,t_high_index,interpolation_factor=1):

        t_fine = np.linspace(self.t[0],self.t[-1],self.t.size*interpolation_factor)
        s0 = interpolate.InterpolatedUnivariateSpline(self.t,u0)
        u0 = s0(t_fine)
        s = interpolate.InterpolatedUnivariateSpline(self.t,u)
        u = s(t_fine)

        window = np.zeros(u.shape,dtype=float)
        window[t_low_index*interpolation_factor:(t_high_index+1)*interpolation_factor] = 1

        u0dot = fftpack.diff(u0,period=self.nt*dt_sec)

        return -(integrate.simps(window*u0dot*(u0-u),dx=dt_sec)/
                integrate.simps(window*u0dot**2,dx=dt_sec))

def fitsread(f):
    return np.squeeze(pyfits.getdata(f)).astype(float)


def find_maxima(arr):
    return np.where(np.diff(np.sign(np.diff(arr)))==-2)[0]



win1 = MainWindow()
win1.connect('destroy', Gtk.main_quit)
win1.show_all()
Gtk.main()
