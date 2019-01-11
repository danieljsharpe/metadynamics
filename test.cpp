/*
Test declaration of class pointer variable with unknown size.
There is a conflict here!
*/

#include <iostream>
using namespace std;

class MyClass {

    public:

    int N;
    int * myptr; // incomplete type (?)
    MyClass(int N);
    void myfunction();
};

MyClass::MyClass(int N1) {

    cout << "In constructor\n";
    N = N1;
    // assign the pointer variable now that its size is known
    int *myptr = new int[N]; // but type is still incomplete
    for (int i=0;i<N;i++) {
        myptr[i] = i;
        cout << myptr[i] << "\n"; }
}

void MyClass::myfunction() {
    cout << "In myfunction\n";
    // but myptr has fallen out of scope
    for (int i=0;i<N;i++) {
        cout << myptr[i] << "\n"; }
}

int main () {

    MyClass myclass1(5);
    myclass1.myfunction();

    return 1.;
}
