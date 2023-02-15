import pygame
import time
from addict import Dict
import traji
from typing import Optional, List, Tuple

INIT_TIME = time.time()

WHITE = (255,255,255)
BLACK = (  0,  0,  0)
RED   = (255,  0,  0)
GREEN = (  0,255,  0)
BLUE  = (  0,  0,255)
YELLOW = (255,255, 0)

# global state storage
STATES = Dict(
    # default: default state with hoving effects
    # drawing: drawing polygon
    stage = None,

    # whether enable extended line when projecting
    extend = False,

    # whether enable projection locally
    local = False,

    # marks events to be handled
    flags = Dict(
        respacing = False,
        densify = False,
    )
)

class Button:
    def __init__(self, text, rect, inactive_color, active_color, action):
        font = pygame.font.Font(None, 30)

        self.rect = pygame.Rect(rect)
        self.text = font.render(text, True, BLACK)
        self.text_rect = self.text.get_rect(center=self.rect.center)
        self.inactive_color = inactive_color
        self.active_color = active_color
        self.is_hover = False
        self.action = action

    def dispatch(self, event):
        if event.type == pygame.MOUSEMOTION:
            self.is_hover = self.rect.collidepoint(event.pos)

        elif event.type == pygame.MOUSEBUTTONDOWN:
            if self.is_hover and self.action:
                self.action()

    def draw(self, screen):
        color = self.active_color if self.is_hover else self.inactive_color
        pygame.draw.rect(screen, color, self.rect)
        screen.blit(self.text, self.text_rect)

def path_to_list(path: traji.Path) -> List[Tuple[float, float]]:
    return [(p.x, p.y) for p in path]

def fib_under(l: float) -> List[int]:
    i, j = 1, 1
    res = [0]
    while j <= l * 1.1:
        i, j = j, i+j
        res.append(i)
    return res

class PolylineDrawer:
    def __init__(self) -> None:
        self.polyline: List[Tuple[float, float]] = []
        self.timestamps: List[float] = [] # for trajectory
        self.last_point: Optional[Tuple[float, float]] = None
        self.last_proj = None
        self.hover_point: Optional[traji.Point] = None
        self.font = pygame.font.Font(None, 20)

    def dispatch(self, event):
        if STATES.stage != 'drawing':
            # handles the changes to the polyline
            if STATES.flags.respacing:
                self.polyline = path_to_list(traji.Path(self.polyline).respacing(10, 50))
                STATES.flags.respacing = False
            if STATES.flags.densify:
                self.polyline = path_to_list(traji.Path(self.polyline).densify(10))
                STATES.flags.densify = False
            if STATES.flags.resample:
                path = traji.Path(self.polyline)
                self.polyline = path_to_list(path.resample_from(fib_under(path.length)))
                STATES.flags.resample = False

            # log the hover position
            if event.type == pygame.MOUSEMOTION:
                self.hover_point = traji.Point(event.pos)
            return

        # initialize
        self.last_proj = None
        if self.last_point is None:
            self.polyline = []
            self.timestamps = []

        if event.type == pygame.MOUSEMOTION:
            self.last_point = event.pos

        elif event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1: # left
                self.polyline.append(event.pos)
                self.timestamps.append(time.time() - INIT_TIME)
            elif event.button == 3: # right
                STATES.stage = 'default'
                self.last_point = None
                print("Exit drawing stage")

    def draw(self, screen):
        if len(self.polyline) == 0:
            return
        elif len(self.polyline) == 1:
            # only has the first point
            if self.last_point is not None:
                pygame.draw.line(screen, WHITE, self.polyline[0], self.last_point)
                pygame.draw.circle(screen, WHITE, self.polyline[0], 3)
        else: # len(self.polyline) > 1
            pygame.draw.lines(screen, WHITE, False, self.polyline + ([self.last_point] if self.last_point else []))
            for p in self.polyline:
                pygame.draw.circle(screen, WHITE, p, 3)

        if STATES.stage != "drawing":
            path = traji.Path(self.polyline)
            if len(self.timestamps) == len(self.polyline):
                # deal with trajectories
                path = traji.Trajectory(path, self.timestamps)
                if STATES.local and (self.last_proj is not None):
                    dist, proj = path.project_local(self.hover_point, self.last_proj, extend=STATES.extend)
                else:
                    dist, proj = path.project(self.hover_point, extend=STATES.extend)
            else:
                # deal with plan paths
                if STATES.local and (self.last_proj is not None):
                    dist, proj = path.project_local(self.hover_point, self.last_proj, extend=STATES.extend)
                else:
                    dist, proj = path.project(self.hover_point, extend=STATES.extend)

            self.last_proj = proj
            perp = path.point_at(proj)
            pygame.draw.line(screen, YELLOW,
                [self.hover_point.x, self.hover_point.y],
                [perp.x, perp.y]
            )

            cx = (self.hover_point.x + perp.x) / 2
            cy = (self.hover_point.y + perp.y) / 2
            text_rect = pygame.Rect((cx - 10, cy - 10, 10, 10))
            if isinstance(path, traji.Trajectory):
                text = self.font.render("%.2f (%.2f) @ %.2f" % (dist, proj.fraction, proj.to_t(path)), True, YELLOW)
            else:
                text = self.font.render("%.2f (%.2f)" % (dist, proj.fraction), True, YELLOW)
            screen.blit(text, text_rect)

class StatusBar:
    def __init__(self) -> None:
        self.font = pygame.font.Font(None, 20)

    def dispatch(self, _event):
        pass # no event will be handled here

    def draw(self, screen):
        text_rect = pygame.Rect((20, 580, 50, 10))
        text = self.font.render("Stage: {}, Extend: {}, Local: {}".format(STATES.stage, STATES.extend, STATES.local),
            True, GREEN)
        screen.blit(text, text_rect)

def get_button_actions(btn):
    if btn == "draw":
        def onclick_draw():
            if STATES.stage == 'drawing':
                STATES.stage = 'default'
                print("Exit drawing stage")
            else:
                STATES.stage = 'drawing'
                print("Enter drawing stage")

        return onclick_draw

    elif btn == "ext":
        def onclick_ext():
            STATES.extend = not STATES.extend
        return onclick_ext

    elif btn == "local":
        def onclick_local():
            STATES.local = not STATES.local
        return onclick_local

    elif btn == "respacing":
        def onclick_respacing():
            STATES.flags.respacing = True
        return onclick_respacing

    elif btn == "densify":
        def onclick_densify():
            STATES.flags.densify = True
        return onclick_densify

    elif btn == "resample":
        def onclick_resample():
            STATES.flags.resample = True
        return onclick_resample

def main():
    
    #### init #####

    pygame.init()
    screen = pygame.display.set_mode((800,600))
    screen_rect = screen.get_rect()

    ##### objects #####

    STATES.stage = "default"
    btn_draw = Button("draw", (20, 20, 80, 30), RED, GREEN, get_button_actions("draw"))
    btn_ext = Button("extend", (20, 60, 80, 30), RED, GREEN, get_button_actions("ext"))
    btn_local = Button("local", (20, 100, 80, 30), RED, GREEN, get_button_actions("local"))
    btn_respacing = Button("respace", (20, 140, 80, 30), RED, GREEN, get_button_actions("respacing"))
    btn_densify = Button("densify", (20, 180, 80, 30), RED, GREEN, get_button_actions("densify"))
    btn_resample = Button("sample", (20, 220, 80, 30), RED, GREEN, get_button_actions("resample"))
    poly_drawer = PolylineDrawer()
    status_bar = StatusBar()
    widgets = [btn_draw, btn_ext, btn_local, btn_respacing, btn_densify, btn_resample, poly_drawer, status_bar]

    ##### mainloop #####

    running = True
    while running:

        ##### events #####
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

            for w in widgets:
                w.dispatch(event)

        ##### draws #####
        screen.fill(BLACK)
        for w in widgets:
            w.draw(screen)
        pygame.display.update()

    #### end ####
    pygame.quit()

if __name__ == "__main__":
    main()