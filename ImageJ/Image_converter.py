from PIL import Image
import os
import sys
import gc
def retrieve_files(path):
    files =[]
    tmp = [filenames for dirpath,dirnames,filenames in os.walk(sys.argv[1])]
    if type(tmp[0]) == type(list()):
        files = [file for file in tmp[0] if 'xls' not in file]
        try:
            os.mkdir('summaries')
        except WindowsError:
            pass
    else:
        files = tmp
    return files
def convert_images(files):
    for file in files:
        gc.collect()
        if "8bit_" in file:
             continue
        im = Image.open(file)
        bw = im.convert("L")
        bw.save('8bit_'+file[:-4]+('.bmp'))
def main():
    files = retrieve_files(sys.argv[1])
    os.chdir(sys.argv[1])
    convert_images(files)
    
main()