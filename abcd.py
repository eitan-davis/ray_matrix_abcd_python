from math import pi
from backend_functions import backend as bd

def lens_thin(f):
    """
    Thin lens 
    
    Parameters
    ----
    f : float
        focal length of lens where f > 0 for convex/positive (converging) lens.
        Only valid if the focal length is much greater than the thickness of the lens. 
    """
    m = [[1,0],[-1/f,1]]
    return bd.array(m)

def lens_thick(n1, n2, r1, r2, t):
    """
    Thick lens
    
    Parameters
    ----
    n1 : float
        refractive index outside of the lens.
    n2 : float
        refractive index of the lens itself (inside the lens).
    R1 : float
        Radius of curvature of First surface.
    R2 : float
        Radius of curvature of Second surface.
    t : float
        center thickness of lens. 
    """
    m_end = bd.array([[1,0],[((n2-n1)/(r2*n1)),(n2/n1)]])
    m_mid = bd.array([[1,t],[0,1]])
    m_start = bd.array([[1,0],[((n1-n2)/(r1*n2)),(n1/n2)]])
    m = m_end @ m_mid @ m_start
    return m

def propagation(d):
    """
    Propagation in free space or in a medium of constant refractive index 

    Parameters
    ----
    d : float
        the distance of the propagation 
    """
    m = [[1,d],[0,1]]
    return bd.array(m)

def mirror():
    """
    Reflection from a flat mirror 
    """
    return bd.eye(2)


def interface(n1, n2):
    """
    Refraction at a flat interface

    Parameters
    ----
    n1 : float
        initial refractive index
    n2 : float
        final refractive index. 

    """
    m = [[1,0],[0,n1/n2]]
    return bd.array(m)

def curved_interface(Re):
    """
    $R_e = R   {cos {\theta}}$ effective radius of curvature in tangential plane (horizontal direction)
    $R_e = R / {cos {\theta}}$ effective radius of curvature in the sagittal plane (vertical direction)
    R = radius of curvature, R > 0 for concave, valid in the paraxial approximation
    {\theta} is the mirror angle of incidence in the horizontal plane. 
    """
    m = [[1,0],[-2/Re,1]]
    return bd.array(m)


def get_new_q(q, m, inverce = False):
    """
    Returns the new q parmeter of the gaussian beam 
    acording to the ABCD matrix and the corrent q parameter.
    The q parameter is in the for of:
    $q = {z - z_0} + i {z_R}$
    Where $z - z_0$ is the relative position to the beam's waist.
    And $z_R$ is the Rayleigh length.

    Parameters
    ----
    q : complex
        The q parameter of the gaussian beam.
    m : 2 * 2 matrix
        the abcd matrix.
    inverce : bool, optional
        it true the inverce value would be calculated
    Returns
    ----
    qt : complex
        The new q parameter of the beam.
    """
    A = m[0][0]
    B = m[0][1]
    C = m[1][0]
    D = m[1][1]
    if not inverce:
        qt = ((A * q + B) / (C * q + D))
    else:
        qt = ((C + D / q) / (A + B / q))
    return qt

def get_gaussian_beams_W0(wavelength :float, z0:float, n = 1.)->float:
  """
  Retuns W0 for a given wave lenght 
  as in Yariv, Amnon (1989). Quantum Electronics (3rd ed.). Wiley. ISBN 0-471-60997-8.

  Parameters:
  ----------
    wavelength : float
      Wave length of the beam
    z0 : float, optional
      Rayleigh range of the gaussian beams. 
    n : float, optional
      Index of refraction, defualt n = 1.
  Retuns:
  ------
    W0 : float
      W0 as in the equetion (the beam's waist radius). 
  """
  return  bd.sqrt(((wavelength / n) * z0 )/ pi)


def get_gaussian_beams_Z0(wavelength :float, w0:float, n = 1.)->float:
  """
  Retuns Rayleigh range for a given wave lenght acording to the beam's waist 
  as in Yariv, Amnon (1989). Quantum Electronics (3rd ed.). Wiley. ISBN 0-471-60997-8.

  Parameters:
  ----------
    wavelength : float
      Wave lenght of the beam
    w0 : float, optional
      The beam's waist radius. 
    n : float, optional
      Index of refraction, defualt n = 1.
  Retuns:
  ------
    Z0 : float
      The Rayleigh range as in the equetion. 
  """
  return  (w0**2 * pi) / ((wavelength / n))