#!/usr/bin/env python3
# -*- mode: Python; coding: utf-8 -*-

"""
@Author: victor
@Date:   19 April 2023 @ 13:18
@Last modified by:   victor
@Last modified by:   victor
@Last modified time: 19 April 2023 @ 13:26

Comment :
"""
import glob, shutil

with open("conf.txt") as f:
    contents = f.readlines()
folder = contents[0].replace("\n", "")

folder += "/*.times"
f.close()
print(f"Folder : {folder}")
for file in glob.glob(folder):
    print(file)
    try:
        shutil.copy(file, file + "x2")
    except:
        pass
