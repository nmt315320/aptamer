function np=DNAcode2(p)
[m,l]=size(p);
for i=1:m
    for j=1:l %j的值是：初始值是1,以2为步长，到py结束
        switch p(i,j) 
            case 0
                np(i,j)= 'C';
            case 1
                np(i,j)= 'T';
            case 2
                np(i,j)= 'A';
            case 3
                np(i,j)= 'G';
        end
    end
end

