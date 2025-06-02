# Kangaroo CPU
## Optimized Pollard's Kangaroo algorithm for ECDSA private key recovery
## Multiprocessing version with significant performance improvements
---

### IMPOSSIBLE. but maybe you're God favorite person.
---
```
pip3 install gmpy2
```

```
$ python3 kangaroo.py 28 03e9e661838a96a65331637e2a3e948dc0756e5009e7cb5c36664d9b72dd18c0a7
[Note] ice_secp256k1 library loaded successfully
[################################################]
[#    Optimized Kangaroo PrivKey Recovery Tool  #]
[#            bitcoin ecdsa secp256k1           #]
[#                  multicore                   #]
[################################################]
[date] Tue Jun  3 02:53:31 2025
[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]
[DEBUG] level=1
[PRNGseed] 1099829438
[library] gmpy2 (optimal performance)
[profile] standart
[range] 2^27..2^28; W = U-L = 0x8000000 (2^28)
[~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~]
[+] Sp-table of pow2 points - ready
[multiprocessing] Running 2 workers...
[✓] Found private key: 0x000000000000000000000000000000000000000000000000000000000d916ce8
[✓] Found private key: 0x000000000000000000000000000000000000000000000000000000000d916ce8

[Results]
[################################################]
Private key found: 0x000000000000000000000000000000000000000000000000000000000d916ce8
Public key: 03e9e661838a96a65331637e2a3e948dc0756e5009e7cb5c36664d9b72dd18c0a7
[################################################]
[date] Tue Jun  3 02:53:32 2025
[exit] done
```
---
##### Donations
- BTC: 3D6VYfj9nCwVXNBcBzTASHdopKGcZzrH1d
  <br><br>
Thank You Very Much
