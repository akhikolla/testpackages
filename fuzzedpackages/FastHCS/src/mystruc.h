struct IdLess {					//internal function.
    template <typename T>
    IdLess(T iter) : values(&*iter) {}
    bool operator()(int left,int right){
        return values[left]<values[right];
    }
    float const* values;
};
