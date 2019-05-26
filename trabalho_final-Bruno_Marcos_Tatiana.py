# -*- coding: utf-8 -*-

# Authors:
#   Bruno Menini Matosak
#   Marcos Antônio de Almeida Rodrigues
#   Tatiana Dias Tardelli Uehara


import gi
gi.require_version('Gtk', '3.0')
from gi.repository import Gtk, GObject
import sys
from matplotlib import pyplot as plt
import wtss
import numpy as np
import datetime
import csv
import scipy.sparse as sparse
from scipy.sparse.linalg import splu
import math


# Global parameters for the functions
# Pyramid filter
window_size_pyramid = 3

# Mean Filter
window_size_mean = 3

# Savitzky-Golay filter
window_size_SG = 5
order_SG = 3
deriv_SG = 0
rate_SG = 1

# Wittaker-Eilers filter
window_size_WT_E = 5
lmbd_WT_E = 10

# Gauss filter
window_size_GA = 5
sigma_GA = 1

# Outliers removal
percent_outliers_removal = 1

class Window(Gtk.ApplicationWindow):

    # Start window
    def __init__(self, app):
        # Iniciates the window
        super(Window, self).__init__(title="MODIS Time Series Filtering", application=app)
        self.set_default_size(500, 700)

        self.grid = Gtk.Grid()

        menubar = Gtk.MenuBar()

        ######
        # Process Menu
        fmi = Gtk.MenuItem.new_with_label("Process")

        menu = Gtk.Menu()
        self.open = Gtk.MenuItem.new_with_label("MOD13Q1")
        self.open.connect("activate", self.data_selection, "MOD13Q1")
        menu.append(self.open)

        self.open = Gtk.MenuItem.new_with_label("MOD13Q1_M")
        self.open.connect("activate", self.data_selection, "MOD13Q1_M")
        menu.append(self.open)

        self.quit = Gtk.MenuItem.new_with_label("Exit")
        self.quit.connect("activate", self.quitApp)
        menu.append(self.quit)

        fmi.set_submenu(menu)

        menubar.add(fmi)

        ######
        # Settings Menu
        fmi = Gtk.MenuItem.new_with_label("Settings")

        menu = Gtk.Menu()

        self.filter_settings = Gtk.MenuItem.new_with_label("Filters")
        self.filter_settings.connect("activate", self.change_filters_settings)
        menu.append(self.filter_settings)

        self.outliers_settings = Gtk.MenuItem.new_with_label("Outliers Removal")
        self.outliers_settings.connect("activate", self.change_outliers_settings)
        menu.append(self.outliers_settings)

        fmi.set_submenu(menu)

        menubar.add(fmi)

        ######
        # Help Menu
        fmi = Gtk.MenuItem.new_with_label("Help")

        menu = Gtk.Menu()

        self.filter_settings = Gtk.MenuItem.new_with_label("About...")
        self.filter_settings.connect("activate", self.about_window)
        menu.append(self.filter_settings)

        self.outliers_settings = Gtk.MenuItem.new_with_label("Help! (still working)")
        self.outliers_settings.connect("activate", self.about_window)
        menu.append(self.outliers_settings)

        fmi.set_submenu(menu)

        menubar.add(fmi)

        self.grid.attach(menubar, 0, 0, 7, 1)
        menubar.set_size_request(500, 5)

        self.add(self.grid)


    # Closes the app when called 'Quit'
    def quitApp(self, par):
        app.quit()

    # Add other widgets to select the parameters to get data from the server
    def data_selection(self, widget, cov):
        self.warning("Note", "Internet connection is required to run this program properly.")
        self.win = Gtk.Window(title="Data: "+cov)
        self.win.connect("destroy", Gtk.main_quit)

        # Latitude entry
        self.entry_lat = Gtk.Entry()
        self.entry_lat.set_text("-22.597063")

        # Longitude entry
        self.entry_long = Gtk.Entry()
        self.entry_long.set_text("-52.221069")

        # Series entry
        list_of_series = Gtk.ListStore(str)
        w = wtss.wtss("http://www.esensing.dpi.inpe.br")
        cv_scheme = w.describe_coverage(cov)
        series = list(cv_scheme["attributes"].keys())
        for serie in series:
            list_of_series.append([serie])
        self.ser_combo = Gtk.ComboBox.new_with_model(list_of_series)
        renderer_text = Gtk.CellRendererText()
        self.ser_combo.pack_start(renderer_text, True)
        self.ser_combo.add_attribute(renderer_text, "text", 0)

        # Start date entry
        self.entry_s_date = Gtk.Entry()
        self.entry_s_date.set_text("2005-01-01")

        # End date entry
        self.entry_e_date = Gtk.Entry()
        self.entry_e_date.set_text("2015-01-01")

        # Pyramid filter check button
        self.check_pyramid = Gtk.CheckButton()
        self.check_pyramid.set_label("Pyramid Filter")

        # Mean filter check button
        self.check_mean = Gtk.CheckButton()
        self.check_mean.set_label("Mean Filter")

        # Gauss filter check button
        self.check_gauss = Gtk.CheckButton()
        self.check_gauss.set_label("Gauss Filter")

        # Savitzky-Golay check button
        self.check_SG = Gtk.CheckButton()
        self.check_SG.set_label("Savitzky-Golay Filter")

        # Wittaker-Eilers check button
        self.check_WT_E = Gtk.CheckButton()
        self.check_WT_E.set_label("Whittaker-Eilers Filter")

        # Filters' Frame
        self.frame_filters = Gtk.Frame()
        self.frame_filters.set_label("Filters")


        hbox = Gtk.Box(orientation=Gtk.Orientation.VERTICAL, spacing=6)
        hbox.pack_start(self.check_pyramid, True, True, 0)
        hbox.pack_start(self.check_mean, True, True, 0)
        hbox.pack_start(self.check_gauss, True, True, 0)
        hbox.pack_start(self.check_SG, True, True, 0)
        hbox.pack_start(self.check_WT_E, True, True, 0)

        self.frame_filters.add(hbox)

        # Info label
        self.info = Gtk.Label()

        # Outliers check button
        self.check_outlier = Gtk.CheckButton()
        self.check_outlier.set_label("Use series with\noutliers removed")

        # Outiers frame
        self.frame_outliers = Gtk.Frame()
        self.frame_outliers.set_label("Outliers")
        self.frame_outliers.add(self.check_outlier)

        # Graph type entry
        graph_type = Gtk.ListStore(str)
        series_g = list(["Line", "Polar"])
        for serie in series_g:
            graph_type.append([serie])
        self.graph_type_combo = Gtk.ComboBox.new_with_model(graph_type)
        renderer_text_graph_type = Gtk.CellRendererText()
        self.graph_type_combo.pack_start(renderer_text_graph_type, True)
        self.graph_type_combo.add_attribute(renderer_text_graph_type, "text", 0)

        # Process button
        self.button = Gtk.Button(label="Process")
        self.button.connect("clicked", self.on_button_clicked, cov, "graph")

        # Process + Save button
        self.button_save = Gtk.Button(label="")
        self.button_save_label = self.button_save.get_child()
        self.button_save_label.set_markup("<small><small>Process + Save to CSV</small></small>")
        self.button_save.connect("clicked", self.on_button_clicked, cov, "save")
        self.button_save.set_property("height-request", 0.5)

        ##### Adicionando os grids #####
        self.grid.attach(Gtk.Label("Lat: "), 0, 1, 1, 1)
        self.grid.attach(self.entry_lat, 1, 1, 1, 1)

        self.grid.attach(Gtk.Label("Long: "), 0, 2, 1, 1)
        self.grid.attach(self.entry_long, 1, 2, 1, 1)

        self.grid.attach(Gtk.Label("Series: "), 0, 3, 1, 1)
        self.grid.attach(self.ser_combo, 1, 3, 1, 1)

        label_s_date = Gtk.Label("Start date:")
        self.grid.attach(label_s_date, 0, 4, 1, 1)
        self.grid.attach(self.entry_s_date, 1, 4, 1, 1)

        self.grid.attach(Gtk.Label("End date:"), 0, 5, 1, 1)
        self.grid.attach(self.entry_e_date, 1, 5, 1, 1)


        self.grid.attach(self.frame_filters, 1, 6, 1, 1)

        self.grid.attach(self.frame_outliers, 1, 7, 1, 1)

        label_s_date = Gtk.Label("Graph Type")
        self.grid.attach(label_s_date, 0, 8, 1, 1)
        self.grid.attach(self.graph_type_combo, 1, 8, 1, 1)


        self.grid.attach(self.button, 1, 10, 1, 1)
        self.grid.attach(self.button_save, 1, 11, 1, 1)
        self.grid.attach_next_to(self.info, self.button, Gtk.PositionType.TOP, 1, 1)
        self.grid.set_row_spacing(6)
        self.grid.set_column_spacing(6)

        self.add(self.grid)

        self.show_all()

    # Things the program do after clicking the button "Process" in the window of the parameters selection
    def on_button_clicked(self, widget, coverage, action):

        tree_iter = self.ser_combo.get_active_iter()
        tree_iter2 = self.graph_type_combo.get_active_iter()
        if tree_iter is not None and tree_iter2 is not None:
            model = self.ser_combo.get_model()
            serie = model[tree_iter][0]

            # Data from the selected point
            try:
                # Check if dates inserted correctly
                date1 = datetime.datetime.strptime(self.entry_s_date.get_text(), "%Y-%m-%d")
                date2 = datetime.datetime.strptime(self.entry_e_date.get_text(), "%Y-%m-%d")
                if (date2-date1).days<0:
                    print((date2-date1).days)
                    int('a')

                # Calls the raw data
                [ti_se, dados] = self.retrieveDataFromPoint(float(self.entry_lat.get_text()),
                                                            float(self.entry_long.get_text()),
                                                            coverage,
                                                            serie,
                                                            self.entry_s_date.get_text(),
                                                            self.entry_e_date.get_text())

                if action == "save":
                    self.get_file( ti_se,
                                   dados,
                                   self.check_pyramid.get_active(),
                                   self.check_mean.get_active(),
                                   self.check_gauss.get_active(),
                                   self.check_outlier.get_active(),
                                   self.check_SG.get_active(),
                                   self.check_WT_E.get_active(),
                                   coverage,
                                   serie,
                                   self.entry_lat.get_text(),
                                   self.entry_long.get_text())
                elif action == "graph":
                    # Calls the funcion to filter if selected and plot the graph
                    self.showGraphFiltered(ti_se,
                                           dados,
                                           self.check_pyramid.get_active(),
                                           self.check_mean.get_active(),
                                           self.check_gauss.get_active(),
                                           self.check_outlier.get_active(),
                                           self.check_SG.get_active(),
                                           self.check_WT_E.get_active(),
                                           serie,
                                           self.graph_type_combo.get_model()[tree_iter2][0])



            except ValueError:
                message = "Input data inserted incorrectly."

                try:
                    float(self.entry_lat.get_text())
                except ValueError:
                    message = message + "\n\t*Lat"

                try:
                    float(self.entry_long.get_text())
                except ValueError:
                    message = message + "\n\t*Long"

                try:
                    datetime.datetime.strptime(self.entry_s_date.get_text(), "%Y-%m-%d")
                except ValueError:
                    message = message + "\n\t*Start date"

                try:
                    datetime.datetime.strptime(self.entry_e_date.get_text(), "%Y-%m-%d")
                except ValueError:
                    message = message + "\n\t*End date"

                try:
                    date1 = datetime.datetime.strptime(self.entry_s_date.get_text(), "%Y-%m-%d")
                    date2 = datetime.datetime.strptime(self.entry_e_date.get_text(), "%Y-%m-%d")

                    if (date2-date1).days<0:
                        message = message + "\n\t*End date < Start date"

                except ValueError:
                    None

                self.error_message("Error", message)

        else:
            self.warning("Warning", "Please select a series and a graph type.")

    # Gets the needed data from the server
    def retrieveDataFromPoint(self, lat, long, series, coverage, date1, date2):
        w = wtss.wtss("http://www.esensing.dpi.inpe.br")
        ts = w.time_series(series, coverage, lat, long, start_date=date1, end_date=date2)
        if len(ts[coverage]) < window_size_mean or (len(ts[coverage]) < window_size_pyramid or (len(ts[coverage]) < window_size_GA or (len(ts[coverage]) < window_size_SG or len(ts[coverage]) < window_size_WT_E))):
            self.warning("Error", "Time interval shorter than what the filter's window sizes allow.")
            return None
        return ts.timeline, ts[coverage]

    # Plots the filtered time series in a new window
    def showGraphFiltered(self, tline, data_raw, f_pyramid, f_mean, f_gauss, f_outlier, f_SG, f_WT_E, coverage, graph_type):

        # Prepares data to print in polar
        def polar(tline,data_raw):
            jdates = tline.copy()
            # julian days
            for i in range(len(jdates)):
                jdates[i] = (datetime.datetime.combine(jdates[i], datetime.datetime.min.time())).timetuple().tm_yday
            jdates = np.asarray(jdates).astype(np.float)
            a = jdates.copy()
            o = jdates.copy()
            for i in range(len(jdates)):
                a[i] = float(data_raw[i] * np.cos(2 * np.pi * jdates[i] / 365.0))
                o[i] = float(data_raw[i] * np.sin(2 * np.pi * jdates[i] / 365.0))
            return a,o

        # Iniciates the plot
        if graph_type=="Line":
            fig, ax = plt.subplots(figsize=(13, 5))
            ax.plot(tline, data_raw, ':', color="grey", label="Raw data")
        elif graph_type=="Polar":
            fig, ax = plt.subplots(figsize=(7, 8))
            [a,o] = polar(tline, data_raw)
            ax.plot(a, o, ':', color="grey", label="Raw data")
        fig.canvas.set_window_title('Graph')

        # Size of the convlution mask
        mask_size = 2

        # Plot a line with outlier removed
        if f_outlier:
            data_filtered = self.remove_outliers(data_raw)
            data_raw = data_filtered.copy()
            if graph_type == "Line":
                ax.plot(tline, data_filtered, '-', linewidth="1", color="grey", label="Outliers Removed")
            elif graph_type=="Polar":
                [a, o] = polar(tline, data_raw)
                ax.plot(a, o, '-', linewidth="1", color="black", label="Outliers Removed")

        # Plot if pyramid filter selected
        if f_pyramid:
            data_filtered = self.filter_pyramid(data_raw)
            if graph_type=="Line":
                ax.plot(tline, data_filtered, '-', color="red", label="Pyramid Filter")
            elif graph_type=="Polar":
                [a, o] = polar(tline, data_filtered)
                ax.plot(a, o, '-', color="red", label="Pyramid Filter")

        # plot if mean filter selected
        if f_mean:
            data_filtered = self.filter_mean(data_raw)
            if graph_type=="Line":
                ax.plot(tline, data_filtered, '-', color="blue", label="Mean Filter")
            elif graph_type=="Polar":
                [a, o] = polar(tline, data_filtered)
                ax.plot(a, o, '-', color="blue", label="Mean Filter")

        # plot if gauss filter selected
        if f_gauss:
            data_filtered = self.filter_gauss(data_raw)
            if graph_type == "Line":
                ax.plot(tline, data_filtered, '-', color="orange", label="Gauss Filter")
            elif graph_type == "Polar":
                [a, o] = polar(tline, data_filtered)
                ax.plot(a, o, '-', color="orange", label="Gauss Filter")

        # plot if Savitzky-Golay selected
        if f_SG:
            data_filtered = self.filter_savitzky_golay(data_raw)
            if graph_type=="Line":
                ax.plot(tline, data_filtered, '-', color="green", label="Savitzky-Golay Filter")
            elif graph_type=="Polar":
                [a, o] = polar(tline, data_filtered)
                ax.plot(a, o, '-', color="green", label="Savitzky-Golay Filter")

        # plot if Savitzky-Golay selected
        if f_WT_E:
            data_filtered = self.filter_whittaker_eilers(data_raw)
            if graph_type == "Line":
                ax.plot(tline, data_filtered, '-', color="purple", label="Wittaker-Eilers Filter")
            elif graph_type == "Polar":
                [a, o] = polar(tline, data_filtered)
                ax.plot(a, o, '-', color="purple", label="Whittaker-Eilers Filter")

        plt.title(coverage.upper() + " Time Series", fontweight='bold')
        if graph_type=="Line":
            ax.set_xlabel("Time", fontweight='bold')
            ax.set_ylabel(coverage.upper() + " value", fontweight='bold')

        plt.legend()
        fig.autofmt_xdate()
        plt.grid()

        plt.show()

    # Opens a dialog so the user can choose a file to save the data
    def get_file(self, tline, data_raw, f_pyramid, f_mean, f_gauss, f_outlier, f_SG, f_WT_E, coverage, series, lat, long):

        # create a filechooserdialog to save:
        # the arguments are: title of the window, parent_window, action,
        # (buttons, response)
        save_dialog = Gtk.FileChooserDialog("Pick a file", self,
                                            Gtk.FileChooserAction.SAVE,
                                            (Gtk.STOCK_CANCEL, Gtk.ResponseType.CANCEL,
                                             Gtk.STOCK_SAVE, Gtk.ResponseType.ACCEPT))
        # the dialog will present a confirmation dialog if the user types a file name that
        # already exists
        save_dialog.set_do_overwrite_confirmation(True)
        # dialog always on top of the textview window
        save_dialog.set_modal(True)
        # connect the dialog to the callback function save_response_cb()
        save_dialog.connect("response", self.save_data, tline, data_raw, f_pyramid, f_mean, f_gauss, f_outlier, f_SG, f_WT_E, coverage, series, lat, long)
        # show the dialog
        save_dialog.show()

    # Process and save data
    def save_data(self, dialog, response_id, tline, data_raw, f_pyramid, f_mean, f_gauss, f_outlier, f_SG, f_WT_E, coverage, series, lat, long):

        save_dialog = dialog
        # if response is "ACCEPT" (the button "Save" has been clicked)
        if response_id == Gtk.ResponseType.ACCEPT:
            # self.file is the currently selected file
            file_name = save_dialog.get_filename()

            if file_name[-4:] != ".csv" and file_name[-4:] != ".CSV":
                file_name = file_name+".csv"

            # Prepares data to print in polar
            def polar(tline, data_raw):
                jdates = tline.copy()
                # julian days
                for i in range(len(jdates)):
                    jdates[i] = (datetime.datetime.combine(jdates[i], datetime.datetime.min.time())).timetuple().tm_yday
                jdates = np.asarray(jdates).astype(np.float)
                a = jdates.copy()
                o = jdates.copy()
                for i in range(len(jdates)):
                    a[i] = float(data_raw[i] * np.cos(2 * np.pi * jdates[i] / 365.0))
                    o[i] = float(data_raw[i] * np.sin(2 * np.pi * jdates[i] / 365.0))
                return a, o

            data_wo_outlier = data_raw.copy()

            # Plot a line with outlier removed
            if f_outlier:
                data_wo_outlier = self.remove_outliers(data_wo_outlier)

            # Plot if pyramid filter selected
            if f_pyramid:
                data_pyramid = self.filter_pyramid(data_wo_outlier)

            # plot if mean filter selected
            if f_mean:
                data_filter_mean = self.filter_mean(data_wo_outlier)

            # plot if mean filter selected
            if f_gauss:
                data_filter_gauss = self.filter_gauss(data_wo_outlier)

            # plot if Savitzky-Golay selected
            if f_SG:
                data_SG = self.filter_savitzky_golay(data_wo_outlier)

            if f_WT_E:
                data_WT_E = self.filter_whittaker_eilers(data_wo_outlier)

            # Opens the CSV file, and writes data
            with open(file_name, mode='w') as CSV:

                CSV_writer = csv.writer(CSV, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)

                # Deals first with the labels
                list = []

                list.append("Coverage:")
                list.append(coverage)
                list.append("Series:")
                list.append(series)
                CSV_writer.writerow(list)

                list = []

                list.append("Latitude:")
                list.append(lat)
                list.append("Longitude:")
                list.append(long)
                CSV_writer.writerow(list)

                list = []
                CSV_writer.writerow(list)

                if f_outlier:
                    list = []
                    list.append("Outliers percentage removal [%]:")
                    global percent_outliers_removal
                    list.append(percent_outliers_removal)
                    CSV_writer.writerow(list)
                if f_pyramid:
                    list = []
                    list.append("Pyramid filter parameters")
                    CSV_writer.writerow(list)
                    list = []
                    list.append("")
                    list.append("Window Size:")
                    global window_size_pyramid
                    list.append(window_size_pyramid)
                    CSV_writer.writerow(list)
                if f_mean:
                    list = []
                    list.append("Mean filter parameters")
                    CSV_writer.writerow(list)
                    list = []
                    list.append("")
                    list.append("Window Size:")
                    global window_size_mean
                    list.append(window_size_mean)
                    CSV_writer.writerow(list)
                if f_gauss:
                    list = []
                    list.append("Gauss filter parameters")
                    CSV_writer.writerow(list)
                    list = []
                    list.append("")
                    list.append("Window Size:")
                    global window_size_GA
                    list.append(window_size_mean)
                    list.append("Standard Deviation:")
                    global sigma_GA
                    list.append(sigma_GA)
                    CSV_writer.writerow(list)
                if f_SG:
                    list = []
                    list.append("Savitzky-Golay filter parameters")
                    CSV_writer.writerow(list)
                    list = []
                    list.append("")
                    global window_size_SG
                    global order_SG
                    global deriv_SG
                    global rate_SG
                    list.append("Window Size:")
                    list.append(window_size_SG)
                    list.append("Polynomial Order:")
                    list.append(order_SG)
                    list.append("Derivative Order:")
                    list.append(deriv_SG)
                    list.append("Rate:")
                    list.append(rate_SG)
                    CSV_writer.writerow(list)
                if f_WT_E:
                    list = []
                    list.append("Whittaker-Eilers filter parameters")
                    CSV_writer.writerow(list)
                    list = []
                    list.append("")
                    global window_size_WT_E
                    global lmbd_WT_E
                    list.append("Window Size:")
                    list.append(window_size_WT_E)
                    list.append("Roughness Penalty:")
                    list.append(order_SG)
                    CSV_writer.writerow(list)


                list = []
                CSV_writer.writerow(list)

                list.append("Date")
                list.append("Data Raw")
                if f_outlier:
                    list.append("Data Without Outliers")
                if f_pyramid:
                    list.append("Pyramid Filter")
                if f_mean:
                    list.append("Mean Filter")
                if f_gauss:
                    list.append("Gauss Filter")
                if f_SG:
                    list.append("Savitzky-Golay Filter")
                if f_WT_E:
                    list.append("Wittaker-Eilers Filter")

                CSV_writer.writerow(list)

                # Writes data itself
                for i in range(len(tline)):
                    list = []

                    list.append(tline[i])
                    list.append((data_raw[i]))

                    if f_outlier:
                        list.append(data_wo_outlier[i])
                    if f_pyramid:
                        list.append(data_pyramid[i])
                    if f_mean:
                        list.append(data_filter_mean[i])
                    if f_gauss:
                        list.append(data_filter_gauss[i])
                    if f_SG:
                        list.append(data_SG[i])
                    if f_WT_E:
                        list.append(data_WT_E[i])

                    CSV_writer.writerow(list)
                # destroy the FileChooserDialog
                dialog.destroy()
                self.warning("Save to CSV complete", "Your file were stored in " + file_name)

        # if response is "CANCEL" (the button "Cancel" has been clicked)
        elif response_id == Gtk.ResponseType.CANCEL:
            print("Save to CSV canceled successfully.")
            dialog.destroy()
        else:
            dialog.destroy()

    #########################################################################
    # - Filtering methods
    #########################################################################
    # pyramid filtering method to the time series. Filter only to the test.
    def filter_pyramid(self, ts):
        ts_filtered = ts.copy()
        global window_size_pyramid
        n_p = window_size_pyramid

        for j in range(n_p, len(ts_filtered) - n_p):
            ts_filtered[j] = 0.
            for i in range(-1 * n_p, n_p + 1):
                ts_filtered[j] = ts_filtered[j] + (abs(abs(i) - (n_p + 1)) * ts[j + i] / ((n_p + 1) ** 2))

        return ts_filtered

    # Mean filtering method
    def filter_mean(self, ts):
        ts_filtered = ts.copy()
        global window_size_mean
        n_p = int((window_size_mean-1)/2)

        for j in range(n_p, len(ts_filtered) - n_p):
            ts_filtered[j] = np.mean(ts[(j-n_p):(j+n_p+1)])

        return ts_filtered

    # Gauss filtering method
    def filter_gauss(self, time_series):

        ts_filtered = time_series.copy()
        global window_size_GA
        global sigma_GA
        mi = 0

        n_p = int((window_size_GA - 1) / 2)

        # defining the mask
        mask = np.zeros(window_size_GA)
        for t in range(-n_p, n_p + 1):
            mask[t + n_p] = float(
                (1 / (sigma_GA * math.sqrt(2 * math.pi))) * math.e ** ((-1 / 2) * (((t - mi) ** 2) / sigma_GA ** 2)))

        for j in range(n_p, len(time_series) - n_p):
            a = np.asarray(time_series[j - n_p:j + n_p + 1])
            soma = np.sum(np.asarray(mask) * a)
            ts_filtered[j] = soma / np.sum(mask)
        return ts_filtered

    # savitzky golay filtering method.
    # Source: https://scipy-cookbook.readthedocs.io/items/SavitzkyGolay.html
    def filter_savitzky_golay(self, ts):
        from math import factorial

        global window_size_SG, order_SG, deriv_SG, rate_SG
        window_size = window_size_SG
        order = order_SG
        deriv = deriv_SG
        rate = rate_SG

        y = np.asarray(ts)

        try:
            window_size = np.abs(np.int(window_size))
            order = np.abs(np.int(order))
        except ValueError:
            print("window_size and order have to be of type int")
        if window_size % 2 != 1 or window_size < 1:
            raise TypeError("window_size size must be a positive odd number")
        if window_size < order + 2:
            raise TypeError("window_size is too small for the polynomials order")
        order_range = range(order + 1)
        half_window = (window_size - 1) // 2
        # precompute coefficients
        b = np.mat([[k ** i for i in order_range] for k in range(-half_window, half_window + 1)])
        m = np.linalg.pinv(b).A[deriv] * rate ** deriv * factorial(deriv)
        # pad the signal at the extremes with
        # values taken from the signal itself
        firstvals = y[0] - np.abs(y[1:half_window + 1][::-1] - y[0])
        lastvals = y[-1] + np.abs(y[-half_window - 1:-1][::-1] - y[-1])
        y = np.concatenate((firstvals, y, lastvals))
        return np.convolve(m[::-1], y, mode='valid')

    # Whittaker filtering method.
    # Source: https://github.com/mhvwerts/whittaker-eilers-smoother/blob/master/whittaker_smooth.py
    def filter_whittaker_eilers(self, y):
        global window_size_WT_E
        global lmbd_WT_E

        d = window_size_WT_E//2 + 1
        lmbd = lmbd_WT_E

        def speyediff(N, d, format='csc'):
            assert not (d < 0), "d must be non negative"
            shape = (N - d, N)
            diagonals = np.zeros(2 * d + 1)
            diagonals[d] = 1.
            for i in range(d):
                diff = diagonals[:-1] - diagonals[1:]
                diagonals = diff
            offsets = np.arange(d + 1)
            spmat = sparse.diags(diagonals, offsets, shape, format=format)
            return spmat

        m = len(y)
        E = sparse.eye(m, format='csc')
        D = speyediff(m, d, format='csc')
        coefmat = E + lmbd * D.conj().T.dot(D)
        z = splu(coefmat).solve(np.asarray(y))

        return z

    # Removing outliers
    # Paper: https://www.mdpi.com/2072-4292/5/12/6159
    def remove_outliers(self, ts):
        ts_corrected = ts.copy()
        global percent_outliers_removal
        per = percent_outliers_removal/100

        for j in range(1, len(ts_corrected) - 1):
            if (ts[j]-ts[j-1]<-per*ts[j-1]) and (ts[j]-ts[j+1]<-per*ts[j+1]):
                ts_corrected[j] = (ts[j-1]+ts[j+1])/2

        return ts_corrected


    #########################################################################
    # - Settings
    #########################################################################
    # change the filters settings
    def change_filters_settings(self, par):
        self.win_edit_filter = Gtk.Window(title="Filters Settings")
        grid_edit = Gtk.Grid()

        # Pyramid
        pyramid_frame = Gtk.Frame()
        pyramid_frame.set_label("Pyramid Filter")
        pyramid_window_label = Gtk.Label("Window Size")
        self.pyramid_entry = Gtk.Entry()
        self.pyramid_entry.set_text(str(window_size_pyramid))
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        hbox.pack_start(pyramid_window_label, True, True, 0)
        hbox.pack_start(self.pyramid_entry, True, True, 0)
        pyramid_frame.add(hbox)
        grid_edit.attach(pyramid_frame, 0, 1, 1, 1)

        # Mean
        mean_frame = Gtk.Frame()
        mean_frame.set_label("Mean Filter")
        mean_window_label = Gtk.Label("Window Size")
        self.mean_entry = Gtk.Entry()
        self.mean_entry.set_text(str(window_size_mean))
        hbox = Gtk.Box(orientation=Gtk.Orientation.HORIZONTAL, spacing=6)
        hbox.pack_start(mean_window_label, True, True, 0)
        hbox.pack_start(self.mean_entry, True, True, 0)
        mean_frame.add(hbox)
        grid_edit.attach(mean_frame, 0, 3, 1, 1)

        # Gauss
        GA_frame = Gtk.Frame()
        GA_frame.set_label("Gauss Filter")
        ####
        grid_GA = Gtk.Grid()
        ####
        self.GA_windowsize_entry = Gtk.Entry()
        self.GA_windowsize_entry.set_text(str(window_size_GA))
        GA_windowsize_label = Gtk.Label("Window Size")
        grid_GA.attach(GA_windowsize_label, 0, 1, 1, 1)
        grid_GA.attach(self.GA_windowsize_entry, 1, 1, 1, 1)
        ####
        self.GA_sigma_entry = Gtk.Entry()
        self.GA_sigma_entry.set_text(str(sigma_GA))
        GA_sigma_label = Gtk.Label("   Stand. Deviation   ")
        grid_GA.attach(GA_sigma_label, 0, 2, 1, 1)
        grid_GA.attach(self.GA_sigma_entry, 1, 2, 1, 1)
        ####
        GA_frame.add(grid_GA)
        grid_edit.attach(GA_frame, 0, 5, 1, 1)

        # SG
        SG_frame = Gtk.Frame()
        SG_frame.set_label("Savitzky-Golay Filter")
        ####
        grid_SG = Gtk.Grid()
        ####
        self.SG_windowsize_entry = Gtk.Entry()
        self.SG_windowsize_entry.set_text(str(window_size_SG))
        SG_windowsize_label = Gtk.Label("Window Size")
        grid_SG.attach(SG_windowsize_label, 0, 1, 1, 1)
        grid_SG.attach(self.SG_windowsize_entry, 1, 1, 1, 1)
        ####
        self.SG_order_entry = Gtk.Entry()
        self.SG_order_entry.set_text(str(order_SG))
        SG_order_label = Gtk.Label("  Polynomial Order  ")
        grid_SG.attach(SG_order_label, 0, 2, 1, 1)
        grid_SG.attach(self.SG_order_entry, 1, 2, 1, 1)
        ####
        self.SG_deriv_entry = Gtk.Entry()
        self.SG_deriv_entry.set_text(str(deriv_SG))
        SG_deriv_label = Gtk.Label("Derivative Order")
        grid_SG.attach(SG_deriv_label, 0, 3, 1, 1)
        grid_SG.attach(self.SG_deriv_entry, 1, 3, 1, 1)
        ####
        self.SG_rate_entry = Gtk.Entry()
        self.SG_rate_entry.set_text(str(rate_SG))
        SG_rate_label = Gtk.Label("Rate")
        grid_SG.attach(SG_rate_label, 0, 4, 1, 1)
        grid_SG.attach(self.SG_rate_entry, 1, 4, 1, 1)
        ####
        SG_frame.add(grid_SG)
        grid_edit.attach(SG_frame, 0, 7, 1, 1)

        # Whittaker-Eilers
        WT_E_frame = Gtk.Frame()
        WT_E_frame.set_label("Whittaker-Eilers Filter")
        ####
        grid_WT_E = Gtk.Grid()
        ####
        self.WT_E_windowsize_entry = Gtk.Entry()
        self.WT_E_windowsize_entry.set_text(str(window_size_WT_E))
        WT_E__windowsize_label = Gtk.Label("Window Size")
        grid_WT_E.attach(WT_E__windowsize_label, 0, 1, 1, 1)
        grid_WT_E.attach(self.WT_E_windowsize_entry, 1, 1, 1, 1)
        ####
        self.WT_E_lmbd_entry = Gtk.Entry()
        self.WT_E_lmbd_entry.set_text(str(lmbd_WT_E))
        WT_E_lmbd_label = Gtk.Label("Roughness Penalty")
        grid_WT_E.attach(WT_E_lmbd_label, 0, 2, 1, 1)
        grid_WT_E.attach(self.WT_E_lmbd_entry, 1, 2, 1, 1)
        ####
        WT_E_frame.add(grid_WT_E)
        grid_edit.attach(WT_E_frame, 0, 9, 1, 1)

        # OK button
        self.bt_change_parameters = Gtk.Button(label="Ok")
        grid_edit.attach(self.bt_change_parameters, 1, 11, 1, 1)
        self.bt_change_parameters.connect("clicked", self.save_filter_new_param)
        self.bt_change_parameters.set_property("height-request", 0.5)

        grid_edit.attach_next_to(Gtk.Label(), pyramid_frame, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label(), mean_frame, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label(), GA_frame, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label(), SG_frame, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label(), WT_E_frame, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label(), self.bt_change_parameters, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label("      "), pyramid_frame, Gtk.PositionType.LEFT, 1, 1)
        grid_edit.attach_next_to(Gtk.Label("      "), self.bt_change_parameters, Gtk.PositionType.RIGHT, 1, 1)
        grid_edit.attach_next_to(Gtk.Label("      "), self.bt_change_parameters, Gtk.PositionType.BOTTOM, 1, 1)

        self.win_edit_filter.add(grid_edit)
        self.win_edit_filter.show_all()

    # Change the outliers settings
    def change_outliers_settings(self, par):
        self.win_edit_outlier = Gtk.Window(title="Outliers Settings")
        # self.win_edit_outlier.set_default_size(410, 430)
        grid_edit = Gtk.Grid()

        self.percent_outlier_label = Gtk.Label("Percentage   ")
        self.percent_outlier_entry = Gtk.Entry()
        self.percent_outlier_entry.set_text(str(rate_SG))

        grid_edit.attach(self.percent_outlier_label, 0, 1, 1, 1)
        grid_edit.attach(self.percent_outlier_entry, 1, 1, 1, 1)

        self.bt_change_outlier = Gtk.Button(label="Ok")
        self.bt_change_outlier.connect("clicked", self.save_outliers_new_param)
        self.bt_change_outlier.set_property("height-request", 0.5)

        grid_edit.attach(self.bt_change_outlier, 2, 3, 1, 1)

        grid_edit.attach_next_to(Gtk.Label(), self.percent_outlier_label, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label("    "), self.percent_outlier_label, Gtk.PositionType.LEFT, 1, 1)
        grid_edit.attach_next_to(Gtk.Label(), self.bt_change_outlier, Gtk.PositionType.TOP, 1, 1)
        grid_edit.attach_next_to(Gtk.Label("    "), self.bt_change_outlier, Gtk.PositionType.RIGHT, 1, 1)
        grid_edit.attach_next_to(Gtk.Label("    "), self.bt_change_outlier, Gtk.PositionType.BOTTOM, 1, 1)


        self.win_edit_outlier.add(grid_edit)
        self.win_edit_outlier.show_all()

    # Saves the new parameters for the filter, defined by the user.
    def save_filter_new_param(self, par):
        # Pyramid filter
        global window_size_pyramid

        # Mean Filter
        global window_size_mean

        # Gauss filter
        global window_size_GA
        global sigma_GA

        # Savitzky-Golay filter
        global window_size_SG
        global order_SG
        global deriv_SG
        global rate_SG

        # Whittaker-Eilers
        global window_size_WT_E
        global lmbd_WT_E


        try:
            window_size_pyramid_test = int(self.pyramid_entry.get_text())
            window_size_mean_test = int(self.mean_entry.get_text())
            window_size_GA_test = int(self.GA_windowsize_entry.get_text())
            window_sigma_GA_test = int(self.GA_sigma_entry.get_text())
            window_size_SG_test = int(self.SG_windowsize_entry.get_text())
            order_SG_test = int(self.SG_order_entry.get_text())
            deriv_SG_test = int(self.SG_deriv_entry.get_text())
            rate_SG_test = int(self.SG_rate_entry.get_text())
            window_size_WT_E_test = int(self.WT_E_windowsize_entry.get_text())
            lmbd_WT_E_test = int(self.WT_E_lmbd_entry.get_text())

            if window_size_pyramid_test % 2 != 1 or window_size_pyramid_test<3:
                int("f")

            if window_size_mean_test % 2 != 1 or window_size_mean_test<3:
                int("f")

            if window_size_GA_test % 2 != 1 or window_size_GA_test<3:
                int("f")

            if window_size_SG_test % 2 != 1 or window_size_SG_test<3:
                int("f")

            if window_size_WT_E_test % 2 != 1 or window_size_WT_E_test<3:
                int("f")

            window_size_pyramid = window_size_pyramid_test
            window_size_mean = window_size_mean_test
            window_size_GA = window_size_GA_test
            sigma_GA = window_sigma_GA_test
            window_size_SG = window_size_SG_test
            order_SG = order_SG_test
            deriv_SG = deriv_SG_test
            rate_SG = rate_SG_test
            window_size_WT_E = window_size_WT_E_test
            lmbd_WT_E = lmbd_WT_E_test

            message = "Values Saved:\n\t-Pyramid Filter\n\t\tWindow Size: "+str(window_size_pyramid)
            message = message + "\n\n\t-Mean Filter\n\t\tWindow Size: "+str(window_size_mean)
            message = message + "\n\n\t-Gauss Filter\n\t\tWindow Size: " + str(window_size_GA)
            message = message + "\n\t\tStandard Deviation: " + str(sigma_GA)
            message = message + "\n\n\t-Savitzky-Golay Filter\n\t\tWindow Size: "+str(window_size_SG)
            message = message + "\n\t\tPolynomial Order: "+str(order_SG)
            message = message + "\n\t\tDerivative Order: "+str(deriv_SG)
            message = message + "\n\t\tRate: "+str(rate_SG)
            message = message + "\n\n\t-Whittaker-Eilers\n\t\tWindow Size: " + str(window_size_WT_E)
            message = message + "\n\t\tRoughness Penalty: " + str(lmbd_WT_E)
            message = message + "\n\nClosing the program will return all parameters to their original values."

            self.warning("Atention", message)
            self.win_edit_filter.destroy()

        except ValueError:
            self.error_message("Error", "Please insert the data correctly. All values must be integer numbers, and the mask orders must be odd numbers equal or greater than 3.")

    # Saves the new parameters to the removal of outliers
    def save_outliers_new_param(self, par):
        global percent_outliers_removal

        try:
            percent_test = float(self.percent_outlier_entry.get_text())

            if percent_test <= 0:
                int("f")

            percent_outliers_removal = percent_test

            message = "Values Saved:\n\t-Outliers Removal\n\t\tPercentage: "+str(percent_outliers_removal)
            message = message + "\n\nClosing the program will return all parameters to their original values."

            self.warning("Atention", message)
            self.win_edit_outlier.destroy()
        except ValueError:
            self.error_message("Error", "Please insert the data correctly. The value must be a number.")

    #########################################################################
    # - Help
    #########################################################################
    def about_window(self, par):
        win = Gtk.Window(title="About")
        grid_about = Gtk.Grid()

        script_title = Gtk.Label()
        script_title.set_markup("<b><big>MODIS (MOD13Q1) Time Series\nFiltering</big></b>")
        script_title.set_justify(Gtk.Justification.CENTER)
        grid_about.attach(script_title, 1,2,1,1)
        grid_about.attach_next_to(Gtk.Label("\n\n"), script_title, Gtk.PositionType.TOP, 1, 1)
        grid_about.attach_next_to(Gtk.Label("\nvers. 0.6.42\n"), script_title, Gtk.PositionType.BOTTOM, 1, 1)
        grid_about.attach_next_to(Gtk.Label("\t\t"), script_title, Gtk.PositionType.LEFT, 1, 1)
        grid_about.attach_next_to(Gtk.Label("\t\t"), script_title, Gtk.PositionType.RIGHT, 1, 1)

        authors = Gtk.Label()
        authors.set_markup("<b>Authors</b>")
        authors.set_justify(Gtk.Justification.LEFT)
        grid_about.attach(authors, 1, 4, 1, 1)

        nomes = Gtk.Label("Bruno Menini Matosak\nMarcos Antônio de Almeida Rodrigues\nTatiana Dias Tardelli Uehara")
        nomes.set_justify(Gtk.Justification.CENTER)
        grid_about.attach(nomes, 1, 5, 1, 1)

        grid_about.attach(Gtk.Label("\n\n"), 1, 6, 1, 1)
        extra = Gtk.Label()
        extra.set_markup("<small>INSTITUTO NACIONAL DE\nPESQUISAS ESPACIAIS\n\nDeveloped for the lecture:\nSER-347 Introdução à Programação\npara Sensoriamento Remoto\n\nProfessors:\nGilberto Ribeiro de Queiroz\nThales Sehn Körting\nFabiano Morelli\n\n</small>")
        grid_about.attach(extra, 1, 7, 1, 1)

        win.add(grid_about)
        win.show_all()


    #########################################################################
    # - Info boxes
    #########################################################################
    # Warning popo-up message
    def warning(self, title, message):
        dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.WARNING,Gtk.ButtonsType.OK, title)
        dialog.format_secondary_text(message)
        dialog.run()

        dialog.destroy()

    def error_message(self, title, message):
        dialog = Gtk.MessageDialog(self, 0, Gtk.MessageType.ERROR, Gtk.ButtonsType.CANCEL, "Error")
        dialog.format_secondary_text(message)
        dialog.run()
        dialog.destroy()

class Application(Gtk.Application):

    def __init__(self):
        super(Application, self).__init__()

    def do_activate(self):
        self.win = Window(self)
        self.win.show_all()

    def do_startup(self):
        Gtk.Application.do_startup(self)

app = Application()
app.run(sys.argv)
