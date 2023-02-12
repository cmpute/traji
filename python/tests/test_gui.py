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
    stage = None
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
        self.hover_point: Optional[traji.Point] = None
        self.font = pygame.font.Font(None, 20)

    def dispatch(self, event):
        if STATES['stage'] != 'drawing':
            if event.type == pygame.MOUSEMOTION:
                self.hover_point = traji.Point(event.pos)
            return

        # initialize
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
            dist, proj = path.project(self.hover_point)
            perp = path.point_at(proj)
            pygame.draw.line(screen, YELLOW,
                [self.hover_point.x, self.hover_point.y],
                [perp.x, perp.y]
            )

            cx = (self.hover_point.x + perp.x) / 2
            cy = (self.hover_point.y + perp.y) / 2
            text_rect = pygame.Rect((cx - 10, cy - 10, 10, 10))
            text = self.font.render("%.2f @ %.2f" % (dist, proj.fraction), True, YELLOW)
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

def main():
    
    #### init #####

    pygame.init()
    screen = pygame.display.set_mode((800,600))
    screen_rect = screen.get_rect()

    ##### objects #####

    STATES['stage'] = "default"
    btn_draw = Button("draw", (20, 20, 50, 30), RED, GREEN, get_button_actions("draw"))
    poly_drawer = PolylineDrawer()
    widgets = [btn_draw, poly_drawer]

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