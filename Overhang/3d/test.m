clear;clc;
x = 0.5*ones(100,50,30);
dc = x;
dv = x;
baseplate = 'S';
[ x, dc,dv ] = AM_filter( x,baseplate, dc, dv );
 