#!/bin/bash
export ARTOOLKIT_CONFIG="filesrc location=6.avi ! avidemux name=demux demux.video_00 ! queue ! decodebin ! ffmpegcolorspace ! videoscale ! videorate ! capsfilter caps=video/x-raw-rgb,bpp=24,width=480,height=320,framerate=30/1 ! identity name=artoolkit ! fakesink"

