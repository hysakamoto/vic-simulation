import time
import sys

class toolbar:

    def __init__(self, toolbar_width, progress=0):
        self.toolbar_width = toolbar_width
        self.progress = progress
        self.update(self.progress)

    def update(self, percentage):
        new_progress = int(percentage*self.toolbar_width)

        if self.progress+0.1 >= new_progress:
            return

        self.progress = new_progress

        # setup toolbar
        sys.stdout.write("[%s]" % (" " * self.toolbar_width))
        sys.stdout.flush()
        sys.stdout.write("\b" * (self.toolbar_width+1)) 

        self.progress = int(percentage*self.toolbar_width)

        for i in xrange(self.progress):
            sys.stdout.write("-")
            sys.stdout.flush()

        sys.stdout.write(" " * (self.toolbar_width-self.progress)+"]")
        sys.stdout.write(" {0:d} %".format(self.progress*100/self.toolbar_width))

        sys.stdout.write("\n")


# tb = toolbar(30)

# for i in xrange(101):
#     time.sleep(0.01) # do real work here
#     # update the bar
#     tb.update(i/100.0)

