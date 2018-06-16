# generate_QAM_baseband
generate the baseband IQ signal for QAM with C++

# build example
```c
// Windows, VS 2015 installed
git clone https://github.com/gongfan99/generate_QAM_baseband.git
cd test
build
test
```

# usage
```c
// Create pool with 3 threads
ThreadPool pool(3);

// Create threads in the pool and let them run
pool.start();

// Submit work
// lambda function is used here as example. Any callable type can be used.
pool.submit( [](float data) { process(data); }, 1.234 );

// Submit work
// it ignores the return value
pool.submitFuture( [](float data) -> float { return process(data); }, 1.234 );

// Submit work that returns future
auto fut = pool.submitFuture( [](float data) -> float { return process(data); }, 1.234 );
float result = fut.get();

// Shutdown the pool, releasing all threads
pool.shutdown();
```

More usage cases can be found in test/test.cc

Either submit() or submitFuture() can be used. submit() may be faster because of less wrapping but it only supports funtion that returns void. submitFuture() supports any form of function but may be slower.

You can come up with an even faster version if the function to be submitted has a fixed known form for example int f(int), then std::bind is not needed in the submit() implementation which can make it faster. 

# Reference
It is modified based on [Mtrebi's thread pool](https://github.com/mtrebi/thread-pool)