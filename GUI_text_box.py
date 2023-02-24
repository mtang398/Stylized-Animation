import PySimpleGUI as sg

def windowInput():
    layout = [  [sg.Text('Enter x_pos: '), sg.InputText()],
            [sg.Text('Enter y_pos: '), sg.InputText()],
            [sg.Button('OK'), sg.Button('Load Recent'), sg.Button('Close'), sg.Button('Save')] ]

    # Create the Window
    window = sg.Window('NURBS Curve GUI', layout)
    # Event Loop to process "events" and get the "values" of the inputs
    while True:
        event, values = window.read()
        if event == sg.WIN_CLOSED or event == 'Close': # if user closes window or clicks close
            window.close()
            return 0, 0, 0
        if event == "OK":
            xpos = int(values[0]) 
            ypos = int(values[1])
            window.close()
            return 1, xpos, ypos
        if event == 'Load Recent':
            window.close()
            return 2,0,0
        if event == 'Save':
            window.close()
            return 4,0,0