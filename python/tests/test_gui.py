import pygame
import traji
from typing import Optional, List, Tuple

WHITE = (255,255,255)
BLACK = (  0,  0,  0)
RED   = (255,  0,  0)
GREEN = (  0,255,  0)
BLUE  = (  0,  0,255)
YELLOW = (255,255, 0)

# global state storage
STATES = dict(
    # default: default state with hoving effects
    # drawing: drawing polygon
    stage = None,

    # whether enable extended line when projecting
    extend = False,

    # whether enable projection locally
    local = False,
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

class PolylineDrawer:
    def __init__(self) -> None:
        self.polyline: List[Tuple[float, float]] = []
        self.last_point: Optional[Tuple[float, float]] = None
        self.last_proj = None
        self.hover_point: Optional[traji.Point] = None
        self.font = pygame.font.Font(None, 20)

    def dispatch(self, event):
        if STATES['stage'] != 'drawing':
            if event.type == pygame.MOUSEMOTION:
                self.hover_point = traji.Point(event.pos)
            return

        # initialize
        self.last_proj = None
        if self.last_point is None:
            self.polyline = []

        if event.type == pygame.MOUSEMOTION:
            self.last_point = event.pos

        elif event.type == pygame.MOUSEBUTTONDOWN:
            if event.button == 1: # left
                self.polyline.append(event.pos)
            elif event.button == 3: # right
                STATES['stage'] = 'default'
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

        if STATES['stage'] != "drawing":
            path = traji.Path(self.polyline)
            if STATES['local'] and (self.last_proj is not None):
                dist, proj = path.project_local(self.hover_point, self.last_proj, extend=STATES['extend'])
            else:
                dist, proj = path.project(self.hover_point, extend=STATES['extend'])
            perp = path.point_at(proj)
            self.last_proj = proj
            pygame.draw.line(screen, YELLOW,
                [self.hover_point.x, self.hover_point.y],
                [perp.x, perp.y]
            )

            cx = (self.hover_point.x + perp.x) / 2
            cy = (self.hover_point.y + perp.y) / 2
            text_rect = pygame.Rect((cx - 10, cy - 10, 10, 10))
            text = self.font.render("%.2f @ %.2f" % (dist, proj.fraction), True, YELLOW)
            screen.blit(text, text_rect)

class StatusBar:
    def __init__(self) -> None:
        self.font = pygame.font.Font(None, 20)

    def dispatch(self, _event):
        pass # no event will be handled here

    def draw(self, screen):
        text_rect = pygame.Rect((20, 580, 50, 10))
        text = self.font.render("Stage: {}, Extend: {}, Local: {}".format(STATES['stage'], STATES['extend'], STATES['local']),
            True, GREEN)
        screen.blit(text, text_rect)

def get_button_actions(btn):
    if btn == "draw":
        def onclick_draw():
            if STATES['stage'] == 'drawing':
                STATES['stage'] = 'default'
                print("Exit drawing stage")
            else:
                STATES['stage'] = 'drawing'
                print("Enter drawing stage")

        return onclick_draw
    elif btn == "ext":
        def onclick_ext():
            STATES['extend'] = not STATES['extend']

        return onclick_ext
    elif btn == "local":
        def onclick_local():
            STATES['local'] = not STATES['local']

        return onclick_local

def main():
    
    #### init #####

    pygame.init()
    screen = pygame.display.set_mode((800,600))
    screen_rect = screen.get_rect()

    ##### objects #####

    STATES['stage'] = "default"
    btn_draw = Button("draw", (20, 20, 80, 30), RED, GREEN, get_button_actions("draw"))
    btn_ext = Button("extend", (20, 60, 80, 30), RED, GREEN, get_button_actions("ext"))
    btn_local = Button("local", (20, 100, 80, 30), RED, GREEN, get_button_actions("local"))
    poly_drawer = PolylineDrawer()
    status_bar = StatusBar()
    widgets = [btn_draw, btn_ext, btn_local, poly_drawer, status_bar]

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