Notes:

EXPLAINING BUILD SCENE:

Build scene is split up into multiple components: 
1) Original scene from A2 (1-37)
   Showcases: Area Light Sources
   How to recreate:
   1) Uncomment this block of code (only)
   2) Run: g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
      Then run: ./RayTracer 512 100 0 areaLS.ppm
   
2) Simple Scene in the water (38-77)
   Showcases: Texture Mapping + Alpha mapping (implemented but not functional)
   How to recreate:
   1) Uncomment this block of code (only)
   2) Run: g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
      Then run: ./RayTracer 512 100 0 textureMap.ppm

3) Simple Dog Scene: Balloon animal dog (78-171)
   Showcases: Antialiasing
   How to recreate:
   1) Uncomment this block of code (only)
   2) Run: g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
      Then run: ./RayTracer 512 100 1 dog_aliasing.ppm
   3) Compare to dog.ppm to see the difference

4) Our Cool Scene: Hello Kitty in her room (172-338)
   Showcases: A cool scene! Hello kitty with more of an abstract murakami theme
   How to recreate:
   1) Uncomment this block of code (only)
   2) Run: g++ -O4 -g svdDynamic.c RayTracer.c utils.c -lm -o RayTracer
      Then run: ./RayTracer 512 100 1 coolKitty.ppm
   




 
