



gcc -I"C:/Program Files/R/R-3.0.1/include" -O4 -Wall  -std=gnu99 -c Gassianintegrals.c -o Gassianintegrals.o
gcc -shared -s -static-libgcc -o Gassianintegrals.dll Gassianintegrals.o -LC:"C:/Program Files/R/R-3.0.1/bin/x64"  

pause