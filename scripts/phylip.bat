@echo off
cd "%2"
if exist outfile del outfile
if exist outtree del outtree
echo Y | "%1" > screenout
del screenout
