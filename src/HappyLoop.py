import HappyHB
import argparse
import os
from os.path import exists
import shutil

def IterateOverFiles(args):

    ext = ('.pdb')
    for file in os.listdir(args.directory):
        if file.endswith(ext):
            print(file)
            AreYouHappy = HappyHB.main(args.directory + '/' + file)
            if AreYouHappy == 0:
                if not exists('Happy'):
                    os.mkdir('Happy')

                shutil.copyfile(args.directory + '/' + file, 'Happy/' + file)
        else:
            continue

parser = argparse.ArgumentParser(description='HappyLoop')
parser.add_argument('--directory', type=str, help='an integer for the accumulator')
args = parser.parse_args()

IterateOverFiles(args)
