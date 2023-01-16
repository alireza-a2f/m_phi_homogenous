class vector {
public:
    double x;
    double y;

    vector();

    void operator=(vector v);
};

vector::vector() {
    x = 0;
    y = 0;
}

void vector::operator=(vector v) {
    this->x = v.x;
    this->y = v.y;
}
