# Build
Install python matplotlib depency (used for ploting)
```
pip3 install matplotlib
```
Build
```
cargo b --release
```

# Run
```
./target/release/calculation problem1 --problem_type=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=u --load_function=type1
```