function colorvec=mycolors(n,col)

green=importdata('green.txt');
blue=importdata('blue.txt');
red=importdata('red.txt');
black=importdata('black.txt');

switch col
    case 'green'
        base=green;
    case 'blue'
        base=blue;
    case 'red'
        base=red;
    case 'black'
        base=black;
end

colorvec=zeros(n,3);
v=0:(1/(n-1)):1;

for i=1:3
    colorvec(:,i)=interp1(0:.01:1,base(:,i),v);
end
end