def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [1, 2, 4, 8, 16, 64, 128, 256 ]
    # print(mem_avail[attempt-1] * 1000, attempt, mem_avail)
    return mem_avail[attempt-1] * 1000