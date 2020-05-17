import tkinter as tk
from tkinter import ttk
class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.master = master
        self.pack()
        self.create_widgets()

    def create_widgets(self):
        self.hi_there = tk.Button(self)
        self.hi_there["text"] = "Hello World\n(click me)"
        self.hi_there["command"] = self.say_hi
        self.hi_there.pack(side="top")
        self.hi_there.pack(expand=0)
        self.contents = tk.StringVar()
        self.contents.set("sdfsdfwefs")
        self.entrythingy = tk.Entry()
        self.entrythingy.pack()
        self.entrythingy["textvariable"] = self.contents
        self.entrythingy.bind('<Key-Return>',
                              self.print_contents)
        self.fred  = tk.Button(self, fg="red", bg="blue",textvariable=self.contents)
        self.fred.pack(side="left")
        self.fred.pack(expand=1)
        self.fred["command"] = self.print_it
        self.fred.bind("<Enter>", self.turn_green)
        self.fred.bind("focus", self.turn_green)
        self.quit = tk.Button(self, text="QUIT", fg="red",
                              command=self.master.destroy)
        self.quit.pack(side="bottom")


    def say_hi(self):
        print("hi there, everyone!")
        # print(self.fred.config(text="11111"))
        self.hi_there["text"] = "222222"
        self.contents.set("qqqqqqqq")
        print(self.contents.get())

    def print_contents(self, event):
        print("hi. contents of entry is now ---->",
              self.contents.get())

    def print_it(self):
        print("hi there")

    def turn_green(self, event):
        event.widget["activeforeground"] = "green"
        self.fred["fg"] = "green"
        print("hi bind")
# tk._test()
# tk.Tcl().eval('info patchlevel')
root = tk.Tk()
app = Application(master=root)
app.mainloop()