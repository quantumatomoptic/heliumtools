import tkinter as tk
import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

def update_plot(*args):
    mean = mean_scale.get()
    std_dev = std_dev_scale.get()
    x = np.linspace(-10, 10, 100)
    y = np.exp(-(x - mean)**2 / (2 * std_dev**2)) / (std_dev * np.sqrt(2 * np.pi))
    
    ax.clear()
    ax.plot(x, y)
    ax.set_xlabel('X')
    ax.set_ylabel('Gaussian Curve')
    ax.set_title('Gaussian Curve Plot')
    
    canvas.draw()

# Create the main window
root = tk.Tk()
root.title('Gaussian Curve Plotter')


# Create sliders for mean and standard deviation
mean_label = tk.Label(root, text='Mean:')
mean_label.pack()
mean_scale = tk.Scale(root, from_=-10, to=10, orient='horizontal', resolution=0.1, command=update_plot)
mean_scale.pack()

std_dev_label = tk.Label(root, text='Standard Deviation:')
std_dev_label.pack()
std_dev_scale = tk.Scale(root, from_=-10, to=10, orient='horizontal', resolution=0.1, command=update_plot)
std_dev_scale.pack()

# Create a Matplotlib figure and canvas to display the plot
fig = Figure()
ax = fig.add_subplot(111)
canvas = FigureCanvasTkAgg(fig, master=root)
canvas.get_tk_widget().pack()

# Start the GUI event loop
root.mainloop()