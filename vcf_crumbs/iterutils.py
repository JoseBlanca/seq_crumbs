
def generate_windows(size, step=None, start=0):
    if step is None:
        step = size

    win_start = None
    while True:
        if win_start is None:
            win_start = start
        else:
            win_start += step
        win_end = win_start + size
        yield win_start, win_end
