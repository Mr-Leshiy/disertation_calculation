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
### Type 1
function u
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=u --problem_type=type1
```
function v
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=v --problem_type=type1
```
function sigma x
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_x --problem_type=type1
```
function sigma y
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_y --problem_type=type1
```
### Type 2
function u
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=u --problem_type=type2
```
function v
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=v --problem_type=type2
```
function sigma x
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_x --problem_type=type2
```
function sigma y
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=sigma_y --problem_type=type2