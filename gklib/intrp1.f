c
c
c
        subroutine intrp1 (a, sizx, sizy, sizz, x, value, errcod)
c ---   Simplest possible interpolation
        implicit none
        integer errcod, sizx, sizy, sizz
        real a(sizx, sizy, sizz), x(3), value

        integer i, j, k
        errcod = 1
        i = nint(x(1))
        j = nint(x(2))
        k = nint(x(3))
        if (i .le. 0) return
        if (j .le. 0) return
        if (k .le. 0) return
        if (i .gt. sizx) return
        if (j .gt. sizy) return
        if (k .gt. sizz) return
        errcod = 0
        value = a(i, j, k)
        return
        end
