import itertools

def ncycle(iterable):
    for item in itertools.cycle(iterable):
        yield item