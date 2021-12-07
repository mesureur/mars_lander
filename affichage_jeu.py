import pygame

#Variables globales
LON = 900
LAR = 500
MA = 100
#Set de couleurs requis
black = (0, 0, 0)

#######         FONCTIONS UTILES         #######
def mise_a_echelle(x, y):
    """Met à l'echelle un point pour qu'il apparaisse au bon endroit
    dans la fenetre de jeu"""
    x = x / 10
    y = y / 10
    #abscisses
    a = LON-MA if x > LON-MA else MA if x < 0 else int(x + MA)
    #ordonnees
    b = LAR-MA if y > LAR-MA else LAR-MA if y < 0 else int(LAR-MA - y)
    
    return a, b


#######         CLASSES         #######
class Fenetre:
    """Initialise la fenetre avec pygame
    Creer l'ecran et affiche les contours de la fenetre de jeu"""
    def __init__(self):
        pygame.init() 
        self.scr = pygame.display.set_mode((LON, LAR), pygame.RESIZABLE)
        pygame.display.set_caption('Mars Lander') 

    def affichage(self):
        self.scr.fill((255, 255, 255))
        pygame.draw.line(self.scr, black, (MA, LAR-MA), (LON-MA, LAR-MA), 4)


class Sol:
    """"Defini le sol de Mars en fonction d'une série de points"""
    def __init__(self, l, screen):
        self.points_sol = [(mise_a_echelle(a[0], a[1])) for a in l]
        self.src = screen

    def affichage(self):
        for i in range (len(self.points_sol)-1):
            pygame.draw.line(self.src, black, self.points_sol[i], self.points_sol[i+1], 2)

    def get_distance(self, point):
        return self.shap_sol.distance(point)

class Text:
    def __init__(self, screen):
        self.src = screen
        self.font = pygame.font.SysFont('arial', 18)

    def affichage(self, nb_population):
        self.text = self.font.render("Génération n°"+str(nb_population),True, black)
        self.src.blit(self.text, (100, 420))