        real  Mols(50,7)

        open(1,file="Mols6.dat")
        open(2,file="Mols6x.dat")

        do 100 ig=1,7
          do 50 j=1,10
             read(1,*)Mols((j-1)*5+1,ig),Mols((j-1)*5+2,ig)
     &,Mols((j-1)*5+3,ig),Mols((j-1)*5+4,ig),Mols((j-1)*5+5,ig)
50        continue
100     continue

        k=50
         write(2,*)k
          do 150 i=1,50
              write(2,*)(Mols(i,ig),ig=1,7)
150        continue
        stop
        end
