# Build
Install python dependencies
```
pip3 install matplotlib
pip3 install argparse
```
Build
```
cargo b --release
```

# Run
## Problem 1
function u
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=u
```
function v
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=v
```
function sigma x
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_x
```
function sigma y
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_y
```
## Problem 2
function u
```
./target/release/calculation problem2 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=u
```
function v
```
./target/release/calculation problem2 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=v
```
function sigma x
```
./target/release/calculation problem2 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_x
```
function sigma y
```
./target/release/calculation problem2 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_y

## Problem 3
function u
```
./target/release/calculation problem3 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --omega=200 --c1=10 --c2=10 --eps=0001 --function_type=u
```
function v
```
./target/release/calculation problem3 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=v
```
function sigma x
```
./target/release/calculation problem3 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_x
```
function sigma y
```
./target/release/calculation problem3 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_y