from tkinter import *
from tkinter import ttk
# root = Tk()
# ttk.Button(root, text="Hello World").grid()
# root.mainloop()

def calculate(*args):
    try:
        value = float(feet.get())
        meters.set((0.3048 * value * 10000.0 + 0.5)/10000.0)
        print(measureSystem)
    except ValueError:
        pass

root = Tk()
root.title("Feet to Meters")

mainframe = ttk.Frame(root, padding="3 3 12 12")
mainframe ['borderwidth'] = 2
mainframe ['relief'] ='sunken'
mainframe.grid(column=0, row=0, sticky=(N, W, E, S))
root.columnconfigure(0, weight=1)
root.rowconfigure(0, weight=1)

feet = StringVar()
meters = StringVar()

feet_entry = ttk.Entry(mainframe, width=7, textvariable=feet)
feet_entry.grid(column=2, row=1, sticky=(W, E))

ttk.Label(mainframe, textvariable=meters).grid(column=2, row=2, sticky=(W, E))
ttk.Button(mainframe, text="Calculate", command=calculate).grid(column=3, row=3, sticky=W)

ttk.Label(mainframe, text="feet").grid(column=3, row=1, sticky=W)
ttk.Label(mainframe, text="is equivalent to").grid(column=1, row=2, sticky=E)
ttk.Label(mainframe, text="meters").grid(column=3, row=2, sticky=W)

l = ttk.Label(mainframe, text="Starting...")
l.grid(column=1, row=4, sticky=(E, W))
l.grid()
l.bind('<Enter>', lambda e: l.configure(text='Moved mouse inside'))
l.bind('<Leave>', lambda e: l.configure(text='Moved mouse outside'))
l.bind('<1>', lambda e: l.configure(text='Clicked left mouse button'))
l.bind('<Double-1>', lambda e: l.configure(text='Double clicked'))
l.bind('<B3-Motion>', lambda e: l.configure(text='right button drag to %d,%d' % (e.x, e.y)))

# lb = Listbox(mainframe)
img_label = ttk.Label(mainframe,text="fenics")
image = PhotoImage(file='fenics.png',height = 200,width=400)
img_label['image'] = image
img_label.grid(column=1, row=5, sticky=(E, W))

button = ttk.Button(mainframe,text ='Okay',command = calculate)
button.state(['disabled'])            # set the disabled flag, disabling the button
# button.state(['!disabled'])           # clear the disabled flag
# button.instate(['disabled'])          # return true if the button is disabled, else false
# button.instate(['!disabled'])         # return true if the button is not disabled, else false
# button.instate(['!disabled'], calculate)    # execute 'cmd' if the button is not disabled

measureSystem = StringVar()
check = ttk.Checkbutton(mainframe, text='Use Metric',
	    command=calculate, variable=measureSystem,
	    onvalue='metric', offvalue='imperial')
check.instate(['alternate'])

phone = StringVar()
home = ttk.Radiobutton(mainframe, text='Home', variable=phone, value='home')
office = ttk.Radiobutton(mainframe, text='Office', variable=phone, value='office')
cell = ttk.Radiobutton(mainframe, text='Mobile', variable=phone, value='cell')

username = StringVar()
name = ttk.Entry(mainframe, textvariable=username)
print('current value is %s' % name.get())
name.delete(0,'end')          # delete between two indices, 0-based
name.insert(0, 'your name')

countryvar = StringVar()
country = ttk.Combobox(mainframe, textvariable=countryvar)
country.bind('<<ComboboxSelected>>', calculate)
country['values'] = ('USA', 'Canada', 'Australia')


for child in mainframe.winfo_children(): child.grid_configure(padx=5, pady=5)

feet_entry.focus()
root.bind('<Return>', calculate)

root.mainloop()