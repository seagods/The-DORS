       Program SplitMols

       integer n
       double precision x(7)
       open(1,file="Mols5.dat")
       open(2,file="Mols5_H2O.dat")
       open(3,file="Mols_CO2.dat")
       open(4,file="Mols5_O3.dat")
       open(5,file="Mols5_N2O.dat")
       open(6,file="Mols5_CO.dat")
       open(7,file="Mols5_CH4.dat")
       open(8,file="Mols_O2.dat")
       read(1,*)n
       write(2,*)n
       write(3,*)n
       write(4,*)n
       write(5,*)n
       write(6,*)n
       write(7,*)n
       write(8,*)n
       do 100 i=1,n
         read(1,*)(x(j),j=1,7)
         write(2,"(1e15.6)")x(1)
         write(3,"(1e15.6)")x(3)
         write(4,"(1e15.6)")x(2)
         write(5,"(1e15.6)")x(5)
         write(6,"(1e15.6)")x(7)
         write(7,"(1e15.6)")x(6)
         write(8,"(1e15.6)")x(4)
100    continue

       stop
       end
