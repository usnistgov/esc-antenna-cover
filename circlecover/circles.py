import circle
import line


def covers_line(circles,l):
    """
    return true if the set of circles covers the line.
    """
    if len(circles) ==  0:
        return False
    c = circles[0]
    t,lines =  c.collides(l)
    newc = list(circles)
    newc.remove(c)
    if not t:
       return covers_line(newc,l)
    if len(lines) == 0:
        # completely enclosed.
        return True
    if len(lines) == 1:
       left = lines[0]
       return covers_line(newc,left)
    elif len(lines) == 2:
       left = lines[0]
       right = lines[1]
       return covers_line(newc,left) and covers_line(newc,right)


                
		
