#!/usr/bin/env python3
# Optimized Pollard's Kangaroo algorithm for ECDSA private key recovery
# Multiprocessing version with significant performance improvements
# https://github.com/rouze-d

from __future__ import print_function
import os
import sys
import time
import math
import random
import ctypes
from multiprocessing import Process, Queue, cpu_count

# Configuration
pow2bits = 28  # Default bits/suborder/exp key
Ntimeit = 2     # Times for avg runtime
timeit_eachnewprvkey = True  # Generate new privkey each loop?
flag_profile = "standart"   # Pollard settings
prngseed = 0    # 0 for random
flag_debug = 1   # Debug level

# Load ice_secp256k1 library
try:
    lib = ctypes.cdll.LoadLibrary("./ice_secp256k1.so")
    from secp256k1 import ice  # Import the ice module from the library
    flag_ice = True
    print("[Note] ice_secp256k1 library loaded successfully")
except Exception as e:
    flag_ice = False
    print(f"[Note] ice_secp256k1 not found - {str(e)}")

# Try to import optimized libraries
try:
    import gmpy2
    from gmpy2 import mpz, invert
    flag_gmpy2 = True
except ImportError:
    flag_gmpy2 = False
    print("[Note] gmpy2 not found - using slower Python math")

# Secp256k1 curve parameters
A_curve = mpz(0) if flag_gmpy2 else 0
B_curve = mpz(7) if flag_gmpy2 else 7
modulo = mpz(0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F) if flag_gmpy2 else 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F
order = mpz(0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141) if flag_gmpy2 else 0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEBAAEDCE6AF48A03BBFD25E8CD0364141
Gx = mpz(0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798) if flag_gmpy2 else 0x79BE667EF9DCBBAC55A06295CE870B07029BFCDB2DCE28D959F2815B16F81798
Gy = mpz(0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8) if flag_gmpy2 else 0x483ADA7726A3C4655DA4FBFC0E1108A8FD17B448A68554199C47D08FFB10D4B8

class Point:
    __slots__ = ['x', 'y']  # Memory optimization
    
    def __init__(self, x=0, y=0):
        self.x = x
        self.y = y

    def __add__(self, other):
        if self.x == 0 and self.y == 0:
            return other
        if other.x == 0 and other.y == 0:
            return self
            
        if self.x == other.x:
            if self.y != other.y:
                return Point(0, 0)
            return self.double()
            
        p = modulo
        dx = (other.x - self.x) % p
        dy = (other.y - self.y) % p
        inv_dx = invert(dx, p) if flag_gmpy2 else pow(dx, p-2, p)
        slope = (dy * inv_dx) % p
        x3 = (slope * slope - self.x - other.x) % p
        y3 = (slope * (self.x - x3) - self.y) % p
        return Point(x3, y3)

    def double(self):
        if self.y == 0:
            return Point(0, 0)
            
        p = modulo
        slope = (3 * self.x * self.x) * invert(2 * self.y, p) % p if flag_gmpy2 else (3 * self.x * self.x) * pow(2 * self.y, p-2, p) % p
        x3 = (slope * slope - 2 * self.x) % p
        y3 = (slope * (self.x - x3) - self.y) % p
        return Point(x3, y3)

    def multiply(self, k):
        if flag_ice:
            # Use ice_secp256k1 for faster multiplication if available
            try:
                result = ice.scalar_multiplication(k)
                return Point(result.x, result.y)
            except:
                pass
                
        # Fall back to standard implementation
        result = Point(0, 0)
        addend = self
        
        while k:
            if k & 1:
                result += addend
            addend = addend.double()
            k >>= 1
        return result

def modinv(a, n=modulo):
    if flag_ice:
        try:
            return ice.modinv(a, n)
        except:
            pass
    return int(invert(mpz(a), mpz(n))) if flag_gmpy2 else pow(a, n-2, n)

def getX2Y(X, y_parity, p=modulo):
    if flag_ice:
        try:
            y = ice.get_y(X, y_parity)
            return y
        except:
            pass
            
    y_squared = (X*X*X + A_curve*X + B_curve) % p
    y = pow(y_squared, (p+1)//4, p)
    if y % 2 != y_parity:
        y = (-y) % p
    return y

def save2file(path, mode, data):
    with open(path, mode) as fp:
        if isinstance(data, (list, tuple, dict, set)):
            fp.writelines(data)
        else:
            fp.write(data)

def time_format(seconds, components=(1,1,1,1,1,1,0,0)):
    intervals = (
        ('y', 31536000), ('m', 2592000), ('d', 86400),
        ('h', 3600), ('min', 60), ('s', 1),
        ('ms', 0.001), ('mcs', 0.000001)
    )
    result = []
    seconds = max(seconds, 0)
    
    for i, (name, count) in enumerate(intervals):
        if components[i]:
            value = int(seconds // count) if count >= 1 else int(seconds / count)
            if value > 0 or (i == 5 and not any(result)):
                seconds -= value * count if count >= 1 else value * count
                result.append(f"{value:02}{name}")
    
    return ' '.join(result) if result else '0s'

def KANGAROOS(W0p, L, U, W, Wsqrt, M):
    G = Point(Gx, Gy)
    Sp = [G]
    for _ in range(255):
        Sp.append(Sp[-1].double())
    
    pow2Jmax = getPow2Jmax(int(round(Wsqrt/2)))
    sizeJmax = 1 << pow2Jmax
    pow2dp = (pow2W//2)-2
    DPmodule = 1 << pow2dp

    dT, dW = M, 1
    Tp, Wp = G.multiply(dT), W0p + G.multiply(dW)
    
    DTp, DWp = {}, {}
    n_jump = last_jump = 0
    prvkey = False
    t0 = t1 = t2 = t1_info = t2_info = time.time()

    while True:
        # Tame kangaroo
        n_jump += 1
        Xcoord = Tp.x % modulo
        
        pw = Xcoord % pow2Jmax
        nowjumpsize = 1 << pw
        
        if Xcoord % DPmodule == 0:
            DTp[Xcoord] = dT
            if flag_debug > 1:
                print(f"\r[tame] T+W={len(DTp)}+{len(DWp)}={len(DTp)+len(DWp)}; {Xcoord:064x} 0x{dT:x}")
                save2file('tame.txt', 'a', f"{Xcoord:064x} {dT}\n")
            
            compare = set(DTp) & set(DWp)
            if compare:
                dDT = DTp[next(iter(compare))]
                dDW = DWp[next(iter(compare))]
                prvkey = abs(dDT - dDW)
                break
        
        dT += nowjumpsize
        Tp += Sp[pw]

        # Wild kangaroo
        n_jump += 1
        Xcoord = Wp.x % modulo
        
        pw = Xcoord % pow2Jmax
        nowjumpsize = 1 << pw
        
        if Xcoord % DPmodule == 0:
            DWp[Xcoord] = dW
            if flag_debug > 1:
                print(f"\r[wild] T+W={len(DTp)}+{len(DWp)}={len(DTp)+len(DWp)}; {Xcoord:064x} 0x{dW:x}")
                save2file('wild.txt', 'a', f"{Xcoord:064x} {dW}\n")
            
            compare = set(DTp) & set(DWp)
            if compare:
                dDT = DTp[next(iter(compare))]
                dDW = DWp[next(iter(compare))]
                prvkey = abs(dDT - dDW)
                break
        
        dW += nowjumpsize
        Wp += Sp[pw]

        # Progress reporting
        if not n_jump % 5000:
            t2 = t2_info = time.time()
            
            if (flag_debug > 0 and (t2_info-t1_info) > 10) or prvkey:
                print(f"\r[i] DP T+W={len(DTp)}+{len(DWp)}={len(DTp)+len(DWp)}; dp/kgr={(len(DTp)+len(DWp))/2:.1f}", end='')
                t1_info = t2_info
            
            if (t2-t1) > 1 or prvkey:
                elapsed = t2-t0
                rate = (n_jump-last_jump)/max(t2-t1, 0.001)
                progress = (n_jump/(2*Wsqrt))*100
                
                print(f"\r[{time_format(elapsed)}; {rate:.0f} j/s; {n_jump}j {progress:.1f}%; dp/kgr={(len(DTp)+len(DWp))/2:.1f}", end='')
                
                if prvkey:
                    timeleft = 0
                else:
                    timeleft = elapsed*(1-(n_jump/(2*Wsqrt)))/(n_jump/(2*Wsqrt))
                
                print(f"; {time_format(timeleft)}]", end='')
                sys.stdout.flush()
                
                t1 = t2
                last_jump = n_jump
        
        if prvkey:
            break

    return prvkey, n_jump, (time.time()-t0), len(DTp), len(DWp)

def getPow2Jmax(optimalmeanjumpsize):
    sumjumpsize = 0
    for i in range(1, 256):
        sumjumpsize += 1 << (i-1)
        now_mean = int(round(sumjumpsize/i))
        next_mean = int(round((sumjumpsize + (1 << i))/i))
        
        if flag_debug > 1:
            print(f"[meanjumpsize#{i:03d}j] {now_mean}(now) <= {optimalmeanjumpsize}(optimal) <= {next_mean}(next)")
        
        if optimalmeanjumpsize - now_mean <= next_mean - optimalmeanjumpsize:
            if flag_debug > 1:
                print(f"[i] pow2Jmax={i} ({now_mean} nearer to optimal)")
            return i
    return 8

def is_on_curve(X, Y, p=modulo):
    return (Y * Y - X * X * X - A_curve * X - B_curve) % p == 0

def getPubkey(prvkey, compressed=True):
    P = Point(Gx, Gy).multiply(prvkey)
    if compressed:
        prefix = '02' if P.y % 2 == 0 else '03'
        return prefix + f"{P.x:064x}"
    else:
        return '04' + f"{P.x:064x}{P.y:064x}"

def usage(bits=28):
    print(f"[usage] python3 {sys.argv[0]} [bits] [pubkey]")
    print(f"        python3 {sys.argv[0]} {bits} 03e9e661838a96a65331637e2a3e948dc0756e5009e7cb5c36664d9b72dd18c0a7")
    print(f"        python3 {sys.argv[0]} 8000000:fffffff 03e9e661838a96a65331637e2a3e948dc0756e5009e7cb5c36664d9b72dd18c0a7")
    print(f"        python3 {sys.argv[0]} 135 02145d2611c823a396ef6712ce0f712f09b9b4f3135e3e0aa3230fb9b6d08d1e16")
    exit(-1)

def kangaroo_worker(pubkey_hex, L, U, result_q, pow2W):
    try:
        X = int(pubkey_hex[2:], 16)
        Y = getX2Y(X, int(pubkey_hex[:2]) - 2)
        if not is_on_curve(X, Y):
            result_q.put((None, 'Invalid curve point'))
            return

        W0p = Point(X, Y)
        W = U - L
        Wsqrt = int(math.isqrt(W))
        M = L + (W // 2)

        prvkey, runjump, runtime, lenT, lenW = KANGAROOS(W0p, L, U, W, Wsqrt, M)
        result_q.put((prvkey, f"jumps={runjump}, time={runtime:.2f}s"))
    except Exception as e:
        result_q.put((None, f"Worker error: {str(e)}"))

if __name__ == '__main__':
    print("[################################################]")
    print("[#    Optimized Kangaroo PrivKey Recovery Tool  #]")
    print("[#            bitcoin ecdsa secp256k1           #]")
    print("[#                  multicore                   #]")
    print("[################################################]")
    print(f"[date] {time.ctime()}")
    print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

    if flag_debug:
        print(f'[DEBUG] level={flag_debug}')

    if not prngseed:
        prngseed = random.randint(1, 1<<32)
    random.seed(prngseed)
    print(f'[PRNGseed] {prngseed}')

    if flag_gmpy2:
        print('[library] gmpy2 (optimal performance)')
    else:
        print('[library] pure Python (slower, install gmpy2 for better performance)')

    print(f'[profile] {flag_profile}')

    # Argument processing
    if len(sys.argv) > 1 and sys.argv[1] in ('--help', '-h', '/?'):
        usage()

    # Initialize range parameters
    if len(sys.argv) > 1:
        try:
            # Try to parse as bit length
            pow2bits = int(sys.argv[1])
            L = 1 << (pow2bits-1)
            U = 1 << pow2bits
            flag_pow2bits = True
            flag_keyspace = False
        except ValueError:
            # Try to parse as key space range
            try:
                L, U = map(lambda x: int(x, 16), sys.argv[1].split(':'))
                flag_pow2bits = False
                flag_keyspace = True
            except:
                usage()
    else:
        # Default to configured pow2bits
        L = 1 << (pow2bits-1)
        U = 1 << pow2bits
        flag_pow2bits = True
        flag_keyspace = False

    # Get public key from command line
    pubkey = None
    if len(sys.argv) > 2:
        pubkey = sys.argv[2]
        if not (pubkey.startswith('02') or pubkey.startswith('03') or pubkey.startswith('04')):
            print("[error] Invalid public key format (must start with 02, 03, or 04)")
            usage()
    else:
        print("[error] Public key argument is required")
        usage()

    W = U - L
    Wsqrt = int(math.isqrt(W))
    M = L + (W // 2)
    pow2W = W.bit_length()

    if flag_pow2bits:
        print(f'[range] 2^{pow2bits-1}..2^{pow2bits}; W = U-L = 0x{W:x} (2^{pow2W})')
    else:
        print(f'[range] 0x{L:x}..0x{U:x}; W = U-L = 0x{W:x} (~2^{pow2W})')

    if pow2W < 8 or pow2W > 160:
        print(f'[error] W must be 2^8..2^160!')
        usage()
    if pow2W > 55:
        print(f'[warn] W = 2^{pow2W} may take a long time to solve')

    print("[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]")

    # Initialize jump table
    Gp = Point(Gx, Gy)
    Zp = Point(0, 0)
    Sp = [Gp]
    for _ in range(255):
        Sp.append(Sp[-1].double())
    print('[+] Sp-table of pow2 points - ready')

    # Main processing loop
    cpu_cores = min(cpu_count(), 4)
    print(f"[multiprocessing] Running {cpu_cores} workers...")

    processes = []
    result_queue = Queue()

    for i in range(cpu_cores):
        p = Process(target=kangaroo_worker, args=(pubkey, L, U, result_queue, pow2W))
        p.start()
        processes.append(p)

    found_key = None
    for p in processes:
        p.join()

    while not result_queue.empty():
        prvkey, status = result_queue.get()
        if prvkey:
            found_key = prvkey
            print(f"[✓] Found private key: 0x{prvkey:064x}")
            save2file('Privkey.txt', 'a', (f"Private key found: 0x{found_key:064x}","\n---------------\n"))

            # Terminate other processes if one found the key
            for proc in processes:
                if proc.is_alive():
                    proc.terminate()
        else:
            print(f"[✗] Worker message: {status}")

    print("\n[Results]")
    print("[################################################]")
    if found_key:
        print(f"Private key found: 0x{found_key:064x}")
        print(f"Public key: {pubkey}")
    else:
        print("No private key found")

    print("[################################################]")
    print(f'[date] {time.ctime()}')
    print('[exit] done')
    exit(0)
