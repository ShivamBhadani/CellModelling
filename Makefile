.RECIPEPREFIX:=>
all:
>gcc cellUKF.c -g -o a.out -DDEBUGG 

run:
>cgdb ./a.out

clean:
>rm ./a.out