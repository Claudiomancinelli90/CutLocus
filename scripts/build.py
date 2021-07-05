#! /usr/bin/env python3 -B
import glob
import os
os.makedirs('build/xcode', exist_ok=True)
os.chdir('build/xcode')
os.system('cmake ../.. -GXcode -DYOCTO_EMBREE=OFF')
os.system('open yocto_gl.xcodeproj')
