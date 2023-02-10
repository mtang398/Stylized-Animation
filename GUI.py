from geomdl import BSpline
from geomdl import knotvector
from pygame import gfxdraw
import pygame,sys
import GUI_text_box
import csv

pygame.init()

# Configurtion set up
class Config():
    def __init__(self):
        self.displaySize = (500,500)
        self.dark = (0,0,0)
        self.grey = (127,127,127)
        self.white = (255,255,255)
        self.red = (255,0,0)
        self.green = (0,255,0)
        self.blue = (0,0,255)
        self.show_points=True

# NURBS Curve set up
class Spline():
    def __init__(self, degree=2):
        self.degree = degree
        self.crv = BSpline.Curve()
        self.crv.degree = degree
        
    def draw(self, points):
        self.crv.ctrlpts = points
        self.crv.knotvector = knotvector.generate(self.crv.degree, self.crv.ctrlpts_size)
        return self.crv.evalpts
        
    def set_degree(self,degree):
        self.crv.degree = degree
        
    def get_curve_points(self):
        return self.crv.evalpts

# Canvas setting        
class Canvas():
    def __init__(self,screen):
        self.screen = screen
        self.curve = Spline()
        self.ctrl_points=[]
        self.count=0
        self.selected=None
        self.move_point=False
        self.add_mode=0
        
        # rect = (x_position, y_position, length, width)
        self.addButton = (220,470,70,20)
        self.deleteButton = (320,470,70,20)
        self.NewButton = (420,470,70,20)
        self.connectButton = (20, 20, 70, 20)
        self.windowButton = (20, 470, 70, 20)
    
    def add_points(self,pts):
        self.ctrl_points.append(pts)
           
    def region(self,button,x,y):
        bx,by,brx,bry=button
        if bx<x<bx+brx and by<y<by+bry:
            return True
        return False       
     
    def button_render(self):
        if self.add_mode:
            pygame.draw.rect(self.screen,cfg.green,self.addButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Add', True, cfg.red, cfg.green)
            textRect = text.get_rect()
            textRect.center = (255, 480)
            screen.blit(text, textRect)
            
            pygame.draw.rect(self.screen,cfg.red,self.windowButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Window', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (55, 480)
            screen.blit(text, textRect)
            
            pygame.draw.rect(self.screen,cfg.red,self.NewButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Clear', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (455, 480)
            screen.blit(text, textRect)
        else:
            pygame.draw.rect(self.screen,cfg.red,self.addButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Add', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (255, 480)
            screen.blit(text, textRect)
            
            pygame.draw.rect(self.screen,cfg.red,self.windowButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Window', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (55, 480)
            screen.blit(text, textRect)
            
            pygame.draw.rect(self.screen,cfg.red,self.NewButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Clear', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (455, 480)
            screen.blit(text, textRect)
            
        if self.selected!=None:
            pygame.draw.rect(self.screen,cfg.red,self.deleteButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Delete', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (355, 480)
            screen.blit(text, textRect)
            
            pygame.draw.rect(self.screen,cfg.red,self.connectButton)
            font = pygame.font.Font('freesansbold.ttf', 12)
            text = font.render('Connect', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (55, 30)
            screen.blit(text, textRect)
            
    def render(self):
        if cfg.show_points:
            if self.count>=2:
                pygame.draw.lines(self.screen, cfg.grey, 0,self.ctrl_points)
            for i,point in enumerate(self.ctrl_points):
                if i==0:
                    pygame.draw.circle(self.screen,cfg.blue,point,5)
                elif i==self.selected:
                    pygame.draw.circle(self.screen,cfg.green,point,5)
                else:
                    pygame.draw.circle(self.screen,cfg.grey,point,5)
        if self.count>=self.curve.degree+1:
            pygame.draw.lines(self.screen,cfg.white,0,self.curve.draw(self.ctrl_points),width=1)
            pts = self.curve.get_curve_points()
            for x in pts:
                pygame.gfxdraw.pixel(screen, int(x[0]), int(x[1]), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]), int(x[1]+1), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]+1), int(x[1]+1), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]-1), int(x[1]+1), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]+1), int(x[1]), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]-1), int(x[1]), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]), int(x[1]-1), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]+1), int(x[1]-1), cfg.red)
                pygame.gfxdraw.pixel(screen, int(x[0]-1), int(x[1]-1), cfg.red)
            # open the file in the write mode
            with open('C:/Users/mstan/project/Potter Project/curve_points_file.csv', 'w', encoding='UTF8') as f:
                # create the csv writer
                writer = csv.writer(f)
                writer.writerow(['x_pos', 'y_pos'])
                writer.writerows(pts)
        self.button_render()
        
    def select(self,x,y,r=10):
        can_add=True
        if self.selected != None:
            if self.region(self.connectButton,x,y):
                self.add_points(self.ctrl_points[0])
                self.count += 1
                self.move_point=True
                can_add=False
        if self.region(self.addButton,x,y):
            self.add_mode=not self.add_mode
            can_add=False
        elif self.region(self.NewButton,x,y):
            self.screen = screen
            self.curve = Spline()
            self.ctrl_points=[]
            self.count=0
            self.selected=None
            self.move_point=False
            self.add_mode=0
            
        elif self.region(self.windowButton,x,y):
            check, xpos, ypos = GUI_text_box.windowInput()
            if check == 1:
                self.add_points((xpos,ypos))
                self.count += 1
                self.move_point=True
                can_add=False
            if check == 2:
                with open('C:/Users/mstan/project/Potter Project/ctrl_points_file.csv', newline='') as f:
                    reader = csv.reader(f)
                    for row in reader:
                        if row != [] and row != ['x_pos', 'y_pos']:
                            self.add_points((int(row[0]),int(row[1])))
                            self.count += 1
                            self.move_point=True
                            can_add=False
            if check == 4: 
                pts = self.ctrl_points
                # open the file in the write mode
                with open('C:/Users/mstan/project/Potter Project/ctrl_points_file.csv', 'w', encoding='UTF8') as f:
                    # create the csv writer
                    writer = csv.writer(f)
                    writer.writerow(['x_pos', 'y_pos'])
                    writer.writerows(pts)
            
        elif self.region(self.deleteButton,x,y):
            self.ctrl_points.pop(self.selected)
            self.count-=1
            self.selected=None
            can_add=False
            
        elif self.count and cfg.show_points:
            for i,points in enumerate(self.ctrl_points):
                px,py = points
                if px-r<x<px+r and py-r<y<py+r:
                    self.selected=i
                    self.move_point=True
                    can_add=False
                    break
        if (self.add_mode and can_add) and cfg.show_points and self.region(self.connectButton,x,y) == False and self.region(self.windowButton,x,y) == False:
            if self.selected != None:
                xpoint,_=self.ctrl_points[self.selected]
                if xpoint>x:
                    self.ctrl_points.insert(self.selected,(x,y))
                else:
                    self.ctrl_points.insert(self.selected+1,(x,y))
            else:
                self.ctrl_points.append((x,y))
            self.count+=1
                    
    def move(self,xy):
        if self.move_point and cfg.show_points:
            self.ctrl_points.pop(self.selected)
            self.ctrl_points.insert(self.selected,xy)
            return True
        
cfg = Config()
screen = pygame.display.set_mode(cfg.displaySize)
canvas = Canvas(screen)

def update():
        screen.fill(cfg.dark)
        #pygame.draw.rect(canvas.screen,cfg.red, canvas.addButton)
        #font = pygame.font.Font('freesansbold.ttf', 12)
        #text = font.render('Add', True, (0,255,0), cfg.red)
        #textRect = text.get_rect()
        #textRect.center = (55, 280)
        #screen.blit(text, textRect)
        #pygame.draw.rect(canvas.screen,cfg.red, canvas.deleteButton)
        #pygame.draw.rect(canvas.screen,cfg.red, canvas.selectButton)
        canvas.render()
        pygame.display.update()
update()

while True:
    for event in pygame.event.get():
            if event.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if event.type == pygame.MOUSEBUTTONDOWN:
                canvas.select(*event.pos)
                update()
            elif event.type == pygame.MOUSEBUTTONUP:
                canvas.move_point=False
            elif event.type == pygame.MOUSEMOTION:
                if canvas.move(event.pos):
                    update()