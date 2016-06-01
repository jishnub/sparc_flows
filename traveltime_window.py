from __future__ import division,print_function
import traveltimes
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

Lx = read_params.get_xlength()
nx = read_params.get_nx()
dt = read_params.get_dt()
dt_min = dt/60
datadir = read_params.get_directory()
codedir = os.path.dirname(os.path.abspath(inspect.getsourcefile(lambda:0)))
masterpixels = np.loadtxt(os.path.join(datadir,"master.pixels"),ndmin=1)
solartime = read_params.get_total_simulation_time()

class MainWindow(Gtk.Window):

    def __init__(self):

        # Files to be used to compute travel times
        self.data_file = None
        self.vzcc_file = None
        self.mode_filter = None
        self.mode_speeds = [0.44,0.60,0.75,0.9,1.2,1.4,1.7,1.9]
        self.source_x = 0
        self.nt = 40+int(solartime*60/dt_min)
        self.t = np.arange(self.nt)*dt
        self.figure = plt.figure()

        Gtk.Window.__init__(self, title="Traveltimes")
        self.set_size_request(400, 150)
        self.set_resizable(False)

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

        self.src_label = Gtk.Label("Source")
        hbox_list[-1].pack_start(self.src_label,expand=True,fill=False,padding=0)

        adjustment_src = Gtk.Adjustment(value=1, lower=1, upper=8,
                                        step_incr=1, page_incr=0, page_size=0)
        self.src = Gtk.SpinButton()
        self.src.set_adjustment(adjustment_src)
        self.src.set_numeric(True)
        self.src.set_snap_to_ticks(True)
        hbox_list[-1].pack_start(self.src,expand=True,fill=False,padding=0)
        self.src.connect("value-changed",self.on_src_changed)

        self.source_x = masterpixels[self.src.get_value_as_int()-1]


        self.iter_label = Gtk.Label("Iter no")
        hbox_list[-1].pack_start(self.iter_label,expand=True,fill=False,padding=0)

        adjustment_iter = Gtk.Adjustment(value=0, lower=0, upper=8,
                                        step_incr=1, page_incr=0, page_size=0)
        self.iter_no = Gtk.SpinButton()
        self.iter_no.set_adjustment(adjustment_iter)
        self.iter_no.set_numeric(True)
        self.iter_no.set_snap_to_ticks(True)
        hbox_list[-1].pack_start(self.iter_no,expand=True,fill=False,padding=0)

        self.get_data_and_vzcc()

        # Load vzcc.fits
        # align_and_hbox()
        #
        # self.vzcc_label = Gtk.Label("vzcc.fits")
        # hbox_list[-1].pack_start(self.vzcc_label,expand=True,fill=False,padding=0)
        #
        # self.vzcc_file_entry = Gtk.Entry()
        # self.vzcc_file_entry.set_width_chars(40)
        # hbox_list[-1].pack_start(self.vzcc_file_entry,expand=True,fill=True,padding=0)
        # self.vzcc_file_entry.connect("changed",self.on_vzcc_entry_changed)
        #
        # self.vzcc_file_browse = Gtk.Button.new_with_label("Browse")
        # hbox_list[-1].pack_start(self.vzcc_file_browse,expand=True,fill=False,padding=0)
        # self.vzcc_file_browse.connect("clicked",self.on_vzcc_file_browse)

        # Load filter
        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        ridge_filters = ["f"]
        ridge_filters.extend(["p"+str(i) for i in xrange(8)])
        self.ridge_filters_list = Gtk.ComboBoxText()
        self.ridge_filters_list.set_entry_text_column(0)
        for ridge in ridge_filters: self.ridge_filters_list.append_text(ridge)
        self.ridge_filters_list.set_active(0)
        self.ridge_filters_list.connect("changed",self.on_filter_change)

        self.filter_label = Gtk.Label("Filter")
        hbox_list[-1].pack_start(self.filter_label,expand=False,fill=False,padding=0)

        hbox_list[-1].pack_start(self.ridge_filters_list,expand=False,fill=False,padding=10)

        self.params_dist = {"min":0,"max":Lx/2}
        self.params_label = Gtk.Label("Min dist: {:.1f}, max dist: {:.1f}".
                                    format(self.params_dist["min"],self.params_dist["max"]))
        self.on_filter_change(self.ridge_filters_list)
        hbox_list[-1].pack_start(self.params_label,expand=False,fill=False,padding=10)

        # Select pixel

        x = np.linspace(-Lx/2,Lx/2,nx,endpoint=False)

        align_and_hbox()
        align_list[-1].set_padding(0,0,10,0)

        self.pixel_label = Gtk.Label("Pixel:")
        hbox_list[-1].pack_start(self.pixel_label,expand=False,fill=False,padding=0)

        adjustment_pix = Gtk.Adjustment(value=0, lower=0, upper=nx-1,
                                        step_incr=1, page_incr=0, page_size=0)
        self.pixel = Gtk.SpinButton()
        self.pixel.set_adjustment(adjustment_pix)
        self.pixel.set_numeric(True)
        self.pixel.set_snap_to_ticks(True)
        self.pixel.connect("value-changed",self.on_pixel_changed)
        hbox_list[-1].pack_start(self.pixel,expand=False,fill=False,padding=0)

        self.x_label = Gtk.Label("x (Mm):")
        hbox_list[-1].pack_start(self.x_label,expand=False,fill=False,padding=0)

        adjustment_x = Gtk.Adjustment(value=x[0], lower=x[0], upper=x[-1],
                                        step_incr=x[1]-x[0], page_incr=0, page_size=0)
        self.xcoord = Gtk.SpinButton()
        self.xcoord.set_adjustment(adjustment_x)
        self.xcoord.set_numeric(True)
        self.xcoord.set_snap_to_ticks(True)
        self.xcoord.connect("value-changed",self.on_xcoord_changed)
        hbox_list[-1].pack_start(self.xcoord,expand=False,fill=False,padding=0)

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

    def on_xcoord_changed(self,xcoord):
        # Check if within params range

        current_x = xcoord.get_value()
        current_pixel = int(nx//2 + current_x/Lx*nx)
        self.pixel.set_value(current_pixel)

    def on_filter_change(self,ridge_filters_list):
        # Try to load params file to get min and max dist
        active_index = ridge_filters_list.get_active()
        params_path = os.path.join(datadir,"params."+str(active_index))
        try:
            params = np.loadtxt(params_path,ndmin=1)
            self.params_dist["min"] = params[0]
            self.params_dist["max"] = params[1]

        except IOError:
            self.params_dist["min"] = 0
            self.params_dist["max"] = Lx/2

        self.params_label.set_text("Min dist: {:.1f}, max dist: {:.1f}".
                                    format(self.params_dist["min"],self.params_dist["max"]))
        filter_file_path = os.path.join(codedir,
                    self.ridge_filters_list.get_active_text()+"mode_filter.fits")
        try:
            self.mode_filter = np.squeeze(pyfits.getdata(filter_file_path)).astype(float)
        except IOError:
            self.mode_filter = None
    #
    # def on_data_file_browse(self,button):
    #     dialog = Gtk.FileChooserDialog("Please choose data.fits", self,
    #         Gtk.FileChooserAction.OPEN,
    #         (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
    #          Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
    #
    #     filter_fits = Gtk.FileFilter()
    #     filter_fits.set_name("FITS files")
    #     filter_fits.add_pattern("*.fits")
    #     dialog.add_filter(filter_fits)
    #
    #     response = dialog.run()
    #     if response == Gtk.ResponseType.OK:
    #         try:
    #             self.data_file = np.squeeze(pyfits.getdata(dialog.get_filename())).astype(float)
    #             self.data_file_entry.set_text(dialog.get_filename())
    #         except IOError:
    #             self.data_file = None
    #             self.data_file_entry.set_text("")
    #             print("Could not load",dialog.get_filename())
    #
    #     dialog.destroy()
    #
    # def on_vzcc_file_browse(self,button):
    #     dialog = Gtk.FileChooserDialog("Please choose vzcc.fits", self,
    #         Gtk.FileChooserAction.OPEN,
    #         (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
    #          Gtk.STOCK_OPEN, Gtk.ResponseType.OK))
    #
    #     filter_fits = Gtk.FileFilter()
    #     filter_fits.set_name("FITS files")
    #     filter_fits.add_pattern("*.fits")
    #     dialog.add_filter(filter_fits)
    #
    #     response = dialog.run()
    #     if response == Gtk.ResponseType.OK:
    #         try:
    #             self.vzcc_file = np.squeeze(pyfits.getdata(dialog.get_filename())).astype(float)
    #             self.vzcc_file_entry.set_text(dialog.get_filename())
    #         except IOError:
    #             self.vzcc_file = None
    #             self.vzcc_file_entry.set_text("")
    #             print("Could not load",dialog.get_filename())
    #
    #     dialog.destroy()

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

    def on_iter_changed(self,widget):
        self.get_data_and_vzcc()

    def get_data_and_vzcc(self):
        # Load data
        try:
            src_file_path = os.path.join(datadir,"tt","data","data"+
                            str(self.src.get_value_as_int()).zfill(2)+".fits")
            self.data_file = fitsread(src_file_path)
            self.nt = self.data_file.shape[0]
            self.time_coordinates = np.arange(self.nt)*dt
        except IOError:
            self.data_file = None

        # Load vzcc
        try:
            vzcc_file_path = os.path.join(datadir,"tt",
                            "iter"+str(self.iter_no.get_value_as_int()).zfill(2),
                            "vzcc_src"+str(self.src.get_value_as_int()).zfill(2)+".fits")
            self.vzcc_file = fitsread(vzcc_file_path)
        except IOError:
            self.vzcc_file = None

    def filter_data_vzcc(self):
        if self.mode_filter is not None:
            data_filtered = np.fft.ifft2(np.fft.fft2(self.data_file)*self.mode_filter)
            vzcc_filtered = np.fft.ifft2(np.fft.fft2(self.vzcc_file)*self.mode_filter)
        else:
            data_filtered = self.data_file
            vzcc_filtered = self.vzcc_file
        return data_filtered,vzcc_filtered

    def on_pixel_tt_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None: return

        pixel = self.pixel.get_value_as_int()

        data_filtered,vzcc_filtered = self.filter_data_vzcc()

        u0 = data_filtered[:,pixel]
        u = vzcc_filtered[:,pixel]
        distance_Mm = abs(self.xcoord.get_value() - self.source_x)
        vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
        timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))
        timefin_index = timest_index + 40

        loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index
        halftime = 10
        t_low_index = loc - halftime
        t_high_index = loc + halftime

        tt_han = compute_tt_hanasoge(u0,u,t_low_index,t_high_index)
        print("tt difference",tt_han,"seconds")
        tt_gb = compute_tt_gizonbirch(u0,u,t_low_index,t_high_index)
        print("tt difference",tt_gb,"seconds")

        plt.cla()
        plt.plot(t/60,u0,label="data",color="blue")
        plt.plot(t/60,u,label="vzcc",color="red")
        plt.axvline(t_low_index*dt_min,ls="dotted",color="black")
        plt.axvline(t_high_index*dt_min,ls="dotted",color="black")
        plt.xlim(timest_index*dt_min,timefin_index*dt_min)
        plt.xlabel("Time (min)")
        plt.ylabel("Wave displacement")
        plt.legend(loc="best")
        plt.show(block=False)

    def on_compute_all_pressed(self,button):
        if self.data_file is None or self.vzcc_file is None: return

        data_filtered,vzcc_filtered = self.filter_data_vzcc()

        tt_hanasoge_list = []
        tt_gb_list = []
        pixel_x_list = []

        for pixel in xrange(nx):

            pixel_x = (pixel/nx-0.5)*Lx
            distance_Mm = abs(pixel_x - self.source_x)
            if distance_Mm<self.params_dist["min"] or distance_Mm>self.params_dist["max"]:
                continue

            vel_Mm_per_min = self.mode_speeds[self.ridge_filters_list.get_active()]
            timest_index = int(np.floor(distance_Mm/vel_Mm_per_min/dt_min))
            timefin_index = timest_index + 40

            u0 = data_filtered[:,pixel]
            u = vzcc_filtered[:,pixel]

            loc = abs(u[timest_index:timefin_index+1]).argmax()+timest_index
            halftime = 10
            t_low_index = loc - halftime
            t_high_index = loc + halftime

            tt = compute_tt_hanasoge(u0,u,t_low_index,t_high_index)
            tt_hanasoge_list.append(tt)
            tt = compute_tt_gizonbirch(u0,u,t_low_index,t_high_index)
            tt_gb_list.append(tt)
            pixel_x_list.append(pixel_x)

        plt.cla()
        plt.plot(pixel_x_list,tt_hanasoge_list,'o',label="Hanasoge")
        plt.plot(pixel_x_list,tt_gb_list,'o',label="Hanasoge")
        plt.legend(loc="best")
        plt.xlabel("x (Mm)")
        plt.ylabel(r"$\Delta \tau$ (sec)")
        plt.show()


def compute_tt_hanasoge(u0, u,t_low_index,t_high_index):
    u0 = u0[t_low_index:t_high_index+1]
    u  =  u[t_low_index:t_high_index+1]
    cc=signal.correlate(u0,u)
    cc_max_index = corr.argmax()
    t = np.arange(-len(u)+1,len(u))*dt

    t_wavepacket = t[cc_max_index-1:cc_max_index+2]
    cc_wavepacket = cc[cc_max_index-1:cc_max_index+2]

    # interpolate to generate more points
    # s = interpolate.UnivariateSpline(t_wavepacket,cc_wavepacket,s=0)
    # t_wavepacket_fine = np.linspace(t_wavepacket[0],t_wavepacket[-1],len(t)*100)
    # cc_wavepacket = s[t_wavepacket]

    p = np.polyfit(t_wavepacket,cc_wavepacket)

    return -p[1]/(2*p[0])

def compute_tt_gizonbirch(u0,u,t_low_index,t_high_index):

    window = np.zeros_like(u)
    window[t_low_index:t_high_index] = 1

    u0dot = fftpack.diff(u0)

    return (integrate.simps(window*u0dot*(u0-u),dx=dt)/
            integrate.simps(window*u0dot**2,dx=dt))

def fitsread(f): return np.squeeze(pyfits.getdata(f)).astype(float)

win1 = MainWindow()
win1.connect("delete-event", Gtk.main_quit)
win1.show_all()
Gtk.main()
