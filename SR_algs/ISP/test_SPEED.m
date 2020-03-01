T1=zeros(1,10000);
T2=zeros(1,10000);

for i=1:10000
    
    tic;y=SCAD_Prox_New(x,1,3);t1=toc;
    tic;y=SCAD_Prox(x,1,3);t2=toc;
    
    T1(i)=t1;
    T2(i)=t2;
    
end