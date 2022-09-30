with open('test.dat', 'w') as f:

    print('hello world', 12, file=f)
    f.write('hello world')

    # f(x) = 3x
    sum = 0
    for i in range(100):
        sum+=1
        print(i)

    print(sum)

    l = [0, 1, 12]
    a = l
    l = []
    a += [34]
    print(a)
    print(l)
