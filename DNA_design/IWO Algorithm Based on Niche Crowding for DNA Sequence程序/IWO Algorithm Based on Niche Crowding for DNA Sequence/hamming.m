function H = hamming(chromlength,x,y)
a = 0;
for j = 1:chromlength
  if x(j) ~= y(j)
     a = a+1;
  end
end
H = a;
    


