#include <iostream>
#include <typeinfo>
#include <cmath>
#include <string>
#include <chrono>
#include <functional>

using namespace std;

template <typename returnType, typename... Args>
struct TimedFunction
{
    function<returnType(Args...)> fn;

    returnType operator()(Args... args)
    {
        chrono::time_point<chrono::system_clock> start, end;
        chrono::duration<double> elapsed;

        start = chrono::system_clock::now();
        std::cout << "running something" << endl;
        returnType result = fn(args...);
        end = chrono::system_clock::now();
        elapsed = end - start;
        printTiming(elapsed);
        return result;
    }

    void printTiming(chrono::duration<double> elapsed)
    {
        const type_info &info = typeid(fn);

        cout << "Executed Call  " << info.name() << ": " << elapsed.count() << " secs" << endl;
    }
};

int add(int x, int y)
{
    return x + y;
}

int main()
{
    using namespace std;

    TimedFunction<int, int, int> timedAdd;
    timedAdd.fn = add;

    int g = timedAdd(4, 5);
    cout << "g = " << g << endl;
    return 0;
}