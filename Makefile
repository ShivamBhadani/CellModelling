.RECIPEPREFIX:=>
all:
>gcc cellUKF.c -lm -g -o a.out

run:
>cgdb ./a.out

clean:
>rm ./a.out