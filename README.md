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
```
./target/release/calculation problem1 --load_function=type1 --a=10 --b=15 --n_x=100 --n_y=10 --puasson_coef=0.25 --young_modulus=200 --eps=0001 --function_type=u --problem_type=type1
```