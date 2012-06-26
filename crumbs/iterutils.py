'''
Created on 2012 eka 26

@author: peio
'''
import random


def sample(iterator, length, sample_size):
    'It makes a sample for the given iterator'

    if not 0 <= sample_size <= length:
        raise ValueError("sample larger than population")

    selected = set(random.randint(0, length - 1) for num in range(sample_size))
    selected_add = selected.add
    while len(selected) < sample_size:
        selected_add(random.randint(0, length - 1))

    for index, item in enumerate(iterator):
        if index in selected:
            yield item


def length(iterator):
    'it counts the items in an iterator. It consumes the iterator'
    count = 0
    for item in iterator:
        count += 1
    return count
