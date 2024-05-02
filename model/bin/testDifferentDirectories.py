import os
import subprocess # just to call an arbitrary command e.g. 'ls'
#ABS_PATH = os.path.dirname(os.path.abspath(__file__))
#print(ABS_PATH)

class cd:
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# Now you can enter the directory like this:
with cd("/root/sequencing/dire/model/dire"):
   # we are in ~/Library
   subprocess.run("ls")

# outside the context manager we are back where we started.
os.getcwd()

