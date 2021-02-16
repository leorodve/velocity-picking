import tkinter
from tkinter import *
from tkinter import ttk
from tkinter import filedialog
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from PIL import ImageTk, Image
from obspy.io.segy.core import _read_segy
from obspy import *
from mpl_toolkits.axes_grid1 import make_axes_locatable
import functions
import colorcet as cc

def upload():
    global filename
    filename = filedialog.askopenfilename(initialdir="/", title="Select A SEG-Y File", filetypes=(("SEG-Y files", "*.segy, *.sgy"),("all files", "*.*")))

def plot_data(plot1, plot2, vm, window, title):    
    figure = Figure(figsize=(10, 4), dpi=100)
    
    ax1 = figure.add_subplot(1, 2, 1)
    ax1.imshow(plot1, extent=[Offsetmin,Offsetmax,(ns-1)*dt,0], cmap='RdBu', vmin=-vm, vmax=vm, aspect='auto')
    ax1.set_xlim(Offsetmin,Offsetmax)
    ax1.set_ylim(ns*dt,0)
    ax1.set_title('Super Gather')
    ax1.set_xlabel('Offset [m]')
    ax1.set_ylabel('Time [ms]')
    
    ax2 = figure.add_subplot(1, 2, 2)
    if title == 'Velocity Spectrum':
        lvl = np.arange(0,max_value, max_value/31)
        b = ax2.contourf(plot2, lvl, origin='upper', extent=[0, 10, ns*dt, 0], cmap=cc.cm.rainbow)
        ax2.set_xticklabels(np.arange(min_vel,max_vel, int((max_vel-min_vel)/5)))
        ax2.set_xlabel('Velocity [m/s]')
        divider = make_axes_locatable(ax2)
        cax = divider.new_horizontal(size="4%", pad=0.03)
        figure.add_axes(cax)
        cbar = figure.colorbar(b, cax=cax)
        cbar.set_ticks([])
    else:
        ax2.imshow(plot2, extent=[Offsetmin,Offsetmax,(ns-1)*dt,0], cmap='RdBu', vmin=-vm, vmax=vm, aspect='auto')
        ax2.set_xlim(Offsetmin,Offsetmax)
        ax2.set_xlabel('Offset [m]')
    ax2.set_ylim(ns*dt,0)
    ax2.set_title(title)
    ax2.set_ylabel('Time [ms]')
    
    canvas = FigureCanvasTkAgg(figure, window)
    return canvas

def offset_calc(gather):
    i = cmp * fold - fold
    offset = np.zeros(fold)
    j = 0
    while j < fold:
        offset[j] = gather.traces[i].stats.segy.trace_header.distance_from_center_of_the_source_point_to_the_center_of_the_receiver_group
        i+=1
        j+=1
    return offset

def mouse_click(event):
    my_label = ttk.Label(root, image=black_square)
    my_label.place(x=event.x-6, y=event.y-6) # place a black square over the position you clicked
    vel_picks.append([event.x, event.y]) #Save the velocity picks to a list

def export_vel_file():
    vel_filename = filedialog.asksaveasfilename(initialdir="/", defaultextension=".*", filetypes=(("Text Files", "*.txt"), ("SEG-Y files", "*.sgy"), ("all files", "*.*")))
    with open(vel_filename,"w") as f:
        f.write("HANDVEL\t%d\n" % cmp)
        f.writelines("%d\t%d\t" % (valueb, valuea) for id7, (valuea, valueb) in enumerate(velocity_picks))

def close_program():
    root.destroy()

def NMO_window():
    global velocity_picks, nmo_mute, nmo_stretch_mute
    for widget in root.winfo_children():
        widget.destroy()
    velocity_function, velocity_picks = functions.create_velocity_function(min_vel, max_vel, ns, dt, vel_picks)
    nmo = functions.nmo_correction(Section, dt_s, Offset, velocity_function)
    nmo_stretch_mute = functions.stretch_mute(velocity_function, ns, dt_s, sm)
    nmo_mute = functions.apply_stretch_mute(nmo_stretch_mute, nmo, fold)
    second_plot = plot_data(Section, nmo, vm, root, 'NMO corrected')
    second_plot.get_tk_widget().grid(row=0, column=0, columnspan=2)
    ttk.Button(root, text="Close program", command=close_program).grid(column=1, row=1, columnspan=1)
    ttk.Button(root, text="Save velocity file", command=export_vel_file).grid(column=0, row=1, columnspan=1)

def welcome_window():
    #labels and grid
    ttk.Button(root, text="Please select a CMP ordered SEG-Y file", command=upload).grid(column=0, row=1, columnspan=3)
    ttk.Label(root, text="Select a CMP").grid(column=0, row=2)
    first_entry = ttk.Entry(root, textvariable=cmp_entry)
    first_entry.grid(column=1, row=2, columnspan=2)
    ttk.Label(root, text="Number of CMP's for super gather").grid(column=0, row=3)
    second_entry = ttk.Entry(root, textvariable=no_cmp)
    second_entry.grid(column=1, row=3, columnspan=2)
    ttk.Label(root, text="Fold").grid(column=0, row=4)
    third_entry = ttk.Entry(root, textvariable=fold_entry)
    third_entry.grid(column=1, row=4, columnspan=2)
    ttk.Label(root, text="Minimum Velocity for CVS").grid(column=0, row=5)
    fourth_entry = ttk.Entry(root, textvariable=min_vel_entry)
    fourth_entry.grid(column=1, row=5, columnspan=2)
    ttk.Label(root, text="Maximum Velocity for CVS").grid(column=0, row=6)
    fifth_entry = ttk.Entry(root, textvariable=max_vel_entry)
    fifth_entry.grid(column=1, row=6, columnspan=2)
    ttk.Label(root, text="Stretch mute").grid(column=0, row=7)
    sixth_entry = ttk.Entry(root, textvariable=stretch_mute_entry)
    sixth_entry.grid(column=1, row=7, columnspan=2)
    ttk.Button(root, text="Next", command=display_data_window).grid(column=0, row=10, columnspan=3)
    
    return root

def display_data_window():
    global ns, dt, dt_s, Section, Offset, vm, Offsetmax, Offsetmin, fold, root, vel, vel_spectrum, cmp, min_vel, max_vel, black_square, vel_picks, max_value, sm
    
    for widget in root.winfo_children():
        widget.destroy()
#   Read segy with obspy
    st = _read_segy(filename)

    nt = len(st.traces)                      #Number of traces in SEGY file
    dt = int(st.traces[0].stats.segy.trace_header.sample_interval_in_ms_for_this_trace/1000) #Sample rate (ms)
    dt_s = dt / 1000                        #Sample rate (s)
    ns = len(st.traces[0].data)             #Number of samples per trace
    cmp = cmp_entry.get()                   #CMP number given by the user
    fold = fold_entry.get()                 #Fold given by the user    
    SG = no_cmp.get() // 2                  #Number of CMP's for super gather
    min_vel = min_vel_entry.get()           #Minimum velocity given by the user
    max_vel = max_vel_entry.get()           #Maximum velocity given by the user
    sm = stretch_mute_entry.get()           #Mute stretch % given by the user

    Section = np.zeros((ns, fold))
    i, k = 0, 0
    while i <= nt-1:
        if st.traces[i].stats.segy.trace_header.ensemble_number == cmp:
            j = int(i // fold - SG)
            a = -SG
            while j <= int((i // fold + SG)) :
                Section[:,k] += st.traces[i + fold * a].data[:]
                j += 1
                a += 1
            i += 1
            k += 1
        else:
            i += 1

    ttk.Label(root, text="CMP: %s"%(cmp)).grid(column=0, row=3)   
    
    vm = np.percentile(Section, 99)
    Offset = offset_calc(st)
    Offsetmax = int(np.amax(Offset))
    Offsetmin = int(np.amin(Offset))
    vel = np.arange(min_vel, max_vel, (max_vel-min_vel)/(11))
    vel_spectrum = functions.calc_velocity_spectrum(Section, dt_s, Offset, vel)    
    max_value = np.amax(vel_spectrum)
    vel_picks = []
    first_plot = plot_data(Section, vel_spectrum, vm, root, 'Velocity Spectrum')
    first_plot.get_tk_widget().grid(row=0, column=0)
    root.bind('<Button-3>', mouse_click)
    black_square = ImageTk.PhotoImage(Image.open("black_square.png"))
    ttk.Button(root, text="Get gather corrected by NMO", command=NMO_window).grid(column=0, row=4, columnspan=3)
    
    return root

#window widget asking the user for file/parameters
root = Tk()
root.title("Velocity Analysis Program")

#variables (default values)
cmp_entry = IntVar()
cmp_entry.set("5")
no_cmp = IntVar()
no_cmp.set("3")
fold_entry = IntVar()
fold_entry.set("10")
min_vel_entry = IntVar()
min_vel_entry.set("800")
max_vel_entry = IntVar()
max_vel_entry.set("5000")
stretch_mute_entry = IntVar()
stretch_mute_entry.set("30")

welcome_frame = welcome_window()
root.mainloop()
