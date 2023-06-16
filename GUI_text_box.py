import PySimpleGUI as sg

# existing functions here

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
import pygame

class InputBox:

    def __init__(self, x, y, w, h, text=''):
        self.rect = pygame.Rect(x, y, w, h)
        self.color = pygame.Color('lightskyblue3')
        self.text = text
        self.txt_surface = pygame.font.Font(None, 32).render(text, True, self.color)
        self.font = pygame.font.Font(None, 32)
        self.active = False

    def handle_event(self, event):
        if event.type == pygame.MOUSEBUTTONDOWN:
            # If the user clicked on the input_box rect.
            if self.rect.collidepoint(event.pos):
                # Toggle the active variable.
                self.active = not self.active
            else:
                self.active = False
        if event.type == pygame.KEYDOWN:
            if self.active:
                if event.key == pygame.K_RETURN:
                    print(self.text)
                    self.text = ''
                elif event.key == pygame.K_BACKSPACE:
                    self.text = self.text[:-1]
                else:
                    self.text += event.unicode

    def update(self, event=None, value=None):
        if event is not None:
            self.text = event
        elif value is not None:
            self.text = value
        else:
            self.text = self.text
        self.txt_surface = self.font.render(self.text, True, self.color)

    def draw(self, screen):
        pygame.draw.rect(screen, self.color, self.rect, 2)
        screen.blit(self.txt_surface, (self.rect.x+5, self.rect.y+5))

    def colorInput(self):
        layout = [  [sg.Text('Enter Red Value: '), sg.InputText()],
                    [sg.Text('Enter Green Value: '), sg.InputText()],
                    [sg.Text('Enter Blue Value: '), sg.InputText()],
                    [sg.Button('OK'), sg.Button('Cancel')] ]

        # Create the Window
        window = sg.Window('RGB Input GUI', layout)
        # Event Loop to process "events" and get the "values" of the inputs
        while True:
            event, values = window.read()
            if event == sg.WIN_CLOSED or event == 'Cancel': # if user closes window or clicks Cancel
                window.close()
                return 0, 0, 0
            if event == "OK":
                r = int(values[0]) 
                g = int(values[1])
                b = int(values[2])
                window.close()
                return r, g, b
