from geomdl import BSpline
from geomdl import knotvector
from pygame import gfxdraw
import pygame,sys
import GUI_text_box
import csv
import os
import time

import numpy as np
from numpy.linalg import solve
import numpy.linalg as linalg

from meshpy.tet import MeshInfo, build
import meshpy.triangle as triangle
import FEM
from scipy.interpolate import LinearNDInterpolator
from matplotlib.path import Path
from numba import njit

pygame.init()

cwd = os.getcwd()

class Config():
    def __init__(self):
        self.displaySize = (500,500)
        self.dark = (0,0,0)
        self.grey = (127,127,127)
        self.white = (255,255,255)
        self.red = (255,0,0)
        self.green = (0,255,0)
        self.blue = (0,0,255)
        self.black = (0, 0, 0)
        self.show_points=True

class Spline():
    def __init__(self, degree=2):
        self.degree = degree
        self.crv = BSpline.Curve()
        self.crv.degree = degree
        
    def draw(self, points):
        self.crv.ctrlpts = [pt['pos'] for pt in points]
        self.crv.knotvector = knotvector.generate(self.crv.degree, self.crv.ctrlpts_size)
        return self.crv.evalpts, self.crv.ctrlpts
        
    def set_degree(self,degree):
        self.crv.degree = degree
        
    def get_curve_points(self):
        return self.crv.evalpts

class Canvas():
    def __init__(self,screen):
        self.screen = screen
        self.curve = Spline()
        self.ctrl_points=[]
        self.count=0
        self.selected=None
        self.move_point=False
        self.add_mode=0
        self.mesh_display = False
        self.mesh = None
        self.fill = False
        self.RGB = []
        self.pts = []
        self.addButton = (220,470,70,20)
        self.deleteButton = (320,470,70,20)
        self.NewButton = (420,470,70,20)
        self.connectButton = (20, 470, 70, 20)
        self.laplace_button = (20, 20, 70, 20)
        self.colorButton = (120,470,70,20)
        
    def get_colors_at_crvpts(self, crvpts):
        crvpts_colors = []
        for crvpt in crvpts:
            x, y = map(int, crvpt)  # convert to integers
            color = self.screen.get_at((x, y))[:-1]  # get the color at the point, omitting the alpha channel

            if color == (0, 0, 0):  # check if the color is the background color
                # Check the surrounding pixels
                for dx, dy in [(-1, 0), (1, 0), (0, -1), (0, 1), (-1, -1), (-1, 1), (1, -1), (1, 1)]:
                    color = self.screen.get_at((x + dx, y + dy))[:-1]
                    if color != (0, 0, 0):  # if a surrounding pixel is not the background color, break the loop
                        break

            crvpts_colors.append(color)
        return crvpts_colors

    def add_points(self,pts,color=(255,255,255)):
        self.ctrl_points.append({'pos':pts,'color':color})
           
    def region(self,button,x,y):
        bx,by,brx,bry=button
        if bx<x<bx+brx and by<y<by+bry:
            return True
        return False
     
    def button_render(self):
        font = pygame.font.Font('freesansbold.ttf', 12)
        if self.add_mode:
            pygame.draw.rect(self.screen,cfg.green,self.addButton)
            text = font.render('Add', True, cfg.red, cfg.green)
            textRect = text.get_rect()
            textRect.center = (255, 480)
            self.screen.blit(text, textRect)

            pygame.draw.rect(self.screen,cfg.red,self.NewButton)
            text = font.render('Clear', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (455, 480)
            self.screen.blit(text, textRect)

        else:
            pygame.draw.rect(self.screen,cfg.red,self.addButton)
            text = font.render('Add', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (255, 480)
            self.screen.blit(text, textRect)

            pygame.draw.rect(self.screen,cfg.red,self.NewButton)
            text = font.render('Clear', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (455, 480)
            self.screen.blit(text, textRect)

        if self.selected!=None:
            pygame.draw.rect(self.screen,cfg.red,self.deleteButton)
            text = font.render('Delete', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (355, 480)
            self.screen.blit(text, textRect)

            pygame.draw.rect(self.screen,cfg.red,self.colorButton)
            text = font.render('Color', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (155, 480)
            self.screen.blit(text, textRect)

            pygame.draw.rect(self.screen,cfg.red,self.connectButton)
            text = font.render('Connect', True, cfg.green, cfg.red)
            textRect = text.get_rect()
            textRect.center = (55, 480)
            self.screen.blit(text, textRect)
            
            pygame.draw.rect(self.screen,cfg.green,self.laplace_button)
            text = font.render('Solve', True, cfg.red, cfg.green)
            textRect = text.get_rect()
            textRect.center = (55, 30)
            self.screen.blit(text, textRect)

        else:
            pygame.draw.rect(self.screen,cfg.green,self.deleteButton)
            text = font.render('Delete', True, cfg.red, cfg.green)
            textRect = text.get_rect()
            textRect.center = (355, 480)
            self.screen.blit(text, textRect)

            pygame.draw.rect(self.screen,cfg.green,self.colorButton)
            text = font.render('Color', True, cfg.red, cfg.green)
            textRect = text.get_rect()
            textRect.center = (155, 480)
            self.screen.blit(text, textRect)

            pygame.draw.rect(self.screen,cfg.green,self.connectButton)
            text = font.render('Connect', True, cfg.red, cfg.green)
            textRect = text.get_rect()
            textRect.center = (55, 480)
            self.screen.blit(text, textRect)
            
    def connect(self):
        self.ctrl_points.append(self.ctrl_points[0].copy())  # Add a copy of the first point to the end
        self.curve.crv.ctrlpts = [pt['pos'] for pt in self.ctrl_points]
        self.curve.crv.knotvector = knotvector.generate(self.curve.crv.degree, len(self.curve.crv.ctrlpts))  # Recalculate the knot vector.
    
    def render(self):
        if self.mesh_display:
            mesh_points = np.array(self.mesh.points)
            mesh_tris = np.array(self.mesh.elements)
            for tri in mesh_tris:
                pygame.draw.polygon(self.screen, cfg.green, [(mesh_points[tri[i]][0], mesh_points[tri[i]][1]) for i in range(3)], width=1)
        if self.fill == True:
            self.color_image(self.pts, self.RGB)
        if cfg.show_points:
            for i,point in enumerate(self.ctrl_points):
                if i==0:
                    pygame.draw.circle(self.screen,cfg.blue,point['pos'],5)
                elif i==self.selected:
                    pygame.draw.circle(self.screen,cfg.green,point['pos'],5)
                else:
                    pygame.draw.circle(self.screen,point['color'],point['pos'],5)
        if self.count >= self.curve.degree+1:
            curve_points = self.curve.draw(self.ctrl_points)[0]
            num_segments = len(self.ctrl_points) - 1

            for i in range(num_segments):
                segment_start = int(i * len(curve_points) / num_segments)
                segment_end = int((i + 1) * len(curve_points) / num_segments)

                for j in range(segment_start, segment_end):
                    t = (j - segment_start) / (segment_end - segment_start)
                    color = (
                        int((1 - t) * self.ctrl_points[i]['color'][0] + t * self.ctrl_points[i + 1]['color'][0]),
                        int((1 - t) * self.ctrl_points[i]['color'][1] + t * self.ctrl_points[i + 1]['color'][1]),
                        int((1 - t) * self.ctrl_points[i]['color'][2] + t * self.ctrl_points[i + 1]['color'][2])
                    )
                    if j < len(curve_points) - 1:
                        pygame.draw.line(self.screen, color, curve_points[j], curve_points[j + 1], width=1)
        self.button_render()
            
    def select(self,x,y,r=10):
        can_add=True
        self.fill = False
        self.RGB = []
        self.pts = []
        if self.selected != None:
            if self.region(self.connectButton,x,y):
                self.connect()
                self.move_point=False
                can_add=False
                self.mesh_display = False
                self.mesh = None
            elif self.region(self.deleteButton,x,y):
                del self.ctrl_points[self.selected]
                self.selected=None
                self.count-=1
                self.move_point=False
                can_add=False
                self.mesh_display = False
                self.mesh = None
            elif self.region(self.colorButton,x,y):
                new_color = GUI_text_box.InputBox.colorInput(self)
                if new_color: # if not cancelled
                    self.ctrl_points[self.selected]['color'] = new_color
                    if self.selected == 0 and len(self.ctrl_points) > 1:  # If the first point's color is changed
                        self.ctrl_points[-1]['color'] = new_color  # Change color of the last point
                    elif self.selected == len(self.ctrl_points) - 1 and len(self.ctrl_points) > 1:  # If the last point's color is changed
                        self.ctrl_points[0]['color'] = new_color  # Change color of the first point
                    can_add=False
            elif self.region(self.laplace_button, x, y):
                self.move_point=False
                can_add=False
                
                print('Solving with FEM ....')
                t1 = time.time()
                curve = self.curve.crv
                crvpts, mesh = self.compute_mesh(curve)
                t2 = time.time()
                print('Time spend creating mesh: ', t2 - t1)
                self.mesh_display = True
                self.mesh = mesh
                
                def find_interpolated_values_inside_curve(curve_points, vertices, values):
                    image_width = 500
                    image_height = 500

                    # Create grid points for interpolation
                    x = np.linspace(0, image_width-1, image_width)
                    y = np.linspace(0, image_height-1, image_height)
                    grid_points = np.meshgrid(x, y)
                    all_points = np.stack(grid_points, axis=-1).reshape(-1, 2)

                    # Check which points are inside the curve
                    polygon = Path(curve_points)
                    inside_curve_mask = polygon.contains_points(all_points)

                    # Select points inside the curve
                    points_inside_curve = all_points[inside_curve_mask]

                    # Create interpolator
                    interpolator = LinearNDInterpolator(vertices, values)

                    # Compute interpolated values for points inside the curve
                    interpolated_values_inside_curve = interpolator(points_inside_curve)

                    return interpolated_values_inside_curve, points_inside_curve

                
                # Calculate and store color values at all the crvpts
                crvpts_colors = self.get_colors_at_crvpts(crvpts)
                R = []
                G = []
                B = []
                
                for r, g, b in crvpts_colors:
                    R.append(r)
                    G.append(g)
                    B.append(b)
                    
                mesh_points = np.array(self.mesh.points)
                mesh_tris = np.array(self.mesh.elements)
                
                RGB_bndry_vals = np.vstack((R, G, B)).T
                
                solutions = FEM.laplace_solver_batched(mesh_points, mesh_tris, crvpts, RGB_bndry_vals)
                t5 = time.time()
                print('Time for solving u: ', t5 - t2)
                
                solutionR, solutionG, solutionB = solutions.T

                pts_inside_crv_R, points_inside_curve = find_interpolated_values_inside_curve(crvpts, mesh_points, solutionR)
                pts_inside_crv_G, points_inside_curve = find_interpolated_values_inside_curve(crvpts, mesh_points, solutionG)
                pts_inside_crv_B, points_inside_curve = find_interpolated_values_inside_curve(crvpts, mesh_points, solutionB)
                t3 = time.time()
                print('Time for three channels of FEM and interpolation: ', t3 - t2)
                
                RGB = []
                for i in range(len(pts_inside_crv_R)):
                    RGB.append((pts_inside_crv_R[i], pts_inside_crv_G[i], pts_inside_crv_B[i]))

                self.fill = True
                self.RGB = RGB
                self.pts = points_inside_curve
                
        if self.count and cfg.show_points:
            for i,points in enumerate(self.ctrl_points):
                px,py = points['pos']
                if (px-r<x<px+r) and (py-r<y<py+r):
                    self.selected = i
                    self.move_point=True
                    can_add=False

        if self.region(self.addButton,x,y):
            self.add_mode = not self.add_mode
            can_add=False
            self.mesh_display = False
            self.mesh = None
        
        if self.region(self.NewButton,x,y):
            self.reset()
            can_add=False
            self.mesh_display = False
            self.mesh = None
        if can_add and self.add_mode:
            self.add_points((x,y))
            self.count+=1
            self.mesh_display = False
            self.mesh = None

    def color_image(self, points, colors):
        surface = pygame.Surface(cfg.displaySize)
        for point, color in zip(points, colors):
            x, y = point
            r, g, b = color
            surface.set_at((int(x), int(y)), (int(r), int(g), int(b)))
        self.screen.blit(surface, (0, 0))

    def compute_mesh(self, curve, N=25, max_area=100):
        def Z(t, curve):
            point = curve.evaluate_single(t)
            return point[:2]  # Discard the last coordinate (if 3D)
        
        def round_trip_connect(start, end):
            return [(i, i + 1) for i in range(start, end)] + [(end, start)]

        def needs_refinement(vertices, area):
            v0, v1, v2 = np.array(vertices)
            dv1, dv2 = v1 - v0, v2 - v0
            n = np.cross(dv1, dv2)
            A = np.linalg.norm(n)/2
            return A > max_area
        
        s_values = np.linspace(0, 1, N+1)
        s_values = s_values[:-1]
        crvpts = [Z(t, curve) for t in s_values]
        
        info = triangle.MeshInfo()
        facets = round_trip_connect(0, len(crvpts) - 1)
        circ_start = len(crvpts)
        facets.extend(round_trip_connect(circ_start, len(crvpts) - 1))
        info.set_points(crvpts)
        info.set_facets(facets)
        mesh = triangle.build(info, refinement_func=needs_refinement)

        return crvpts, mesh

    def move(self,x,y):
        if self.move_point:
            self.ctrl_points[self.selected]['pos']=(x,y)
            if self.mesh != None:
                curve = self.curve.crv
                crvpts, mesh = self.compute_mesh(curve)

                self.mesh_display = True
                self.mesh = mesh

    def reset(self):
        self.ctrl_points=[]
        self.count=0
        self.selected=None
        self.move_point=False
        self.add_mode=0
        self.mesh_display = False
        self.mesh = None

cfg = Config()
screen = pygame.display.set_mode(cfg.displaySize)
canvas = Canvas(screen)
running=True
clock = pygame.time.Clock()

while running:
    clock.tick(60)
    screen.fill(cfg.dark)
    canvas.render()
    pygame.display.flip()

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False
        if event.type == pygame.MOUSEBUTTONDOWN:
            x, y = pygame.mouse.get_pos()
            canvas.select(x,y)
        if event.type == pygame.MOUSEBUTTONUP:
            canvas.move_point = False
        if event.type == pygame.MOUSEMOTION:
            x, y = pygame.mouse.get_pos()
            canvas.move(x,y)
