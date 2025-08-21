import sys
import time
import datetime
import os

export_path = sys.argv[1] if len(sys.argv) > 1 else "resource.txt"
interval = sys.argv[2] if len(sys.argv) > 2 else 1

if os.path.exists(export_path):
    os.remove(export_path)

#write a function that reads processor use % from /proc/stat
def read_proc():
    with open('/proc/stat', 'r') as f:
        lines = f.readlines()
    
    cpu_line = lines[0].strip().split()
    total_time = sum(int(x) for x in cpu_line[1:])  # Sum of all time fields
    idle_time = int(cpu_line[4])  # Idle time is the 5th field

    # Calculate CPU usage percentage
    cpu_usage = (total_time - idle_time) / total_time * 100 if total_time > 0 else 0
    return round(cpu_usage, 2)

#write a function that reads memory use % from /proc/meminfo
def read_mem():
    with open('/proc/meminfo', 'r') as f:
        lines = f.readlines()
    mem_total = 0
    mem_free = 0
    for line in lines:
        if line.startswith('MemTotal:'):
            mem_total = int(line.split()[1])
        elif line.startswith('MemAvailable:'):
            mem_free = int(line.split()[1])
    # Calculate memory usage percentage
    mem_usage = ((mem_total - mem_free) / mem_total * 100) if mem_total > 0 else 0
    return round(mem_usage,2)

#make a loop that runs every interval seconds and writes the cpu and memory usage to the file
def main():
    while True:
        cpu_usage = read_proc()
        mem_usage = read_mem()
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        with open(export_path, "a") as f:
            f.write(f"{timestamp}, {cpu_usage}, {mem_usage}\n")
        time.sleep(1)

if __name__ == "__main__":
    main()
