
from tkinter import *
from tkinter import ttk as ttk



# def _Preview(*args, **kwargs):
#     global CalcType, MethodType

def _CalcType(*args, **kwargs):
    global CalcType
    Calc = CalcType.get()
    print(Calc)

def _MethodType(*args, **kwargs):
    global MethodType
    Method  = MethodType.get()
    print(Method)
    if MethodType.get() == 'DFT':
        global FuncType

        # Creation of Functional tpype dropdown box
        Label(MethodTab, text='Functional', font = ("Times New Roman", 12)).grid(column=3, row=1, padx=10, pady=25)

        FuncType = StringVar()

        FunctionalTypes = ttk.Combobox(MethodTab, textvariable=FuncType, values=('M062-X', '\u03c9B97X-D', 'PBE', 'PBE0'), state='readonly')
        FunctionalTypes.current()
        FunctionalTypes.grid(row=1 ,column=4)
        FunctionalTypes.bind("<<ComboboxSelected>>", _FunctionalType)

        

def _FunctionalType(*args, **kwargs):
    global FuncType
    Functional = FuncType.get()
    print(Functional)


if __name__ == '__main__':
    root = Tk()
    root.title('XYZ to DALTON')
    tabControl = ttk.Notebook(root)

    # Tabs can be added here
    MainTab = Frame(tabControl)
    MethodTab = Frame(tabControl)

    tabControl.add(MainTab, text='Main')
    tabControl.add(MethodTab, text='Method')
    tabControl.pack(expand=1, fill="both")

    # Creation of Job type dropdown box (combobox)
    Label(MainTab, text='Job Type', font = ("Times New Roman", 12)).grid(column=0, row=1, padx=10, pady=25)

    CalcType = StringVar()

    CalculationTypes = ttk.Combobox(MainTab, textvariable=CalcType, values=('Energy', 'Optimization', 'Opt+Freq', 'Frequency', 'IRC'), state='readonly')
    CalculationTypes.current(0)
    CalculationTypes.grid(row=1, column=1)
    CalculationTypes.bind("<<ComboboxSelected>>", _CalcType)
    
    # Creation of Method type dropdown box
    Label(MethodTab, text='Method', font = ("Times New Roman", 12)).grid(column=0, row=1, padx=10, pady=25)

    MethodType = StringVar()

    MethodTypes = ttk.Combobox(MethodTab, textvariable=MethodType, values=('HF', 'MP2', 'CCSD', 'DFT'), state='readonly')
    MethodTypes.current(0)
    MethodTypes.grid(row=1, column=1)
    MethodTypes.bind("<<ComboboxSelected>>", _MethodType)

    root.mainloop()