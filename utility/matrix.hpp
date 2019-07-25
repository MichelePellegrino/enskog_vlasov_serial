#ifndef EV_MASK_MATRIX_HPP
#define EV_MASK_MATRIX_HPP

#include <Eigen/Dense>
#include <cassert>
#include <array>
#include <algorithm>

#define TL 0    /*!< top-left       */
#define CL 1    /*!< centre-left    */
#define BL 2    /*!< bottom-left    */
#define BC 3    /*!< bottom-centre  */
#define BR 4    /*!< bottom-right   */
#define CR 5    /*!< centre-right   */
#define TR 6    /*!< top-right      */
#define TC 7    /*!< top-center     */

#define N_BUF 8   /*!< total number of buffers */

namespace ev_matrix {

  using namespace Eigen;

  // ###################################################################### //
  // ##### MASK MATRIX #################################################### //
  // ###################################################################### //

  template <class data_type, int StorageOrder = RowMajor>
  class MaskMatrix : public Array<data_type, Dynamic, Dynamic, StorageOrder>
  {
  protected:
    typedef Array<data_type, Dynamic, Dynamic, StorageOrder> DynamicMatrix;
    typedef Block<const DynamicMatrix> constDynamicBlock;
    int low_x, up_x, low_y, up_y;
    // Index transformations
    inline int idx_i(const int& i) const
    {
      return i - low_x;
    }
    inline int idx_j(const int& j) const
    {
      return j - low_y;
    }
  public:
    MaskMatrix ( int lx, int ux, int ly, int uy ):
      DynamicMatrix(ux-lx, uy-ly), low_x(lx), up_x(ux), low_y(ly), up_y(uy)
      {
        assert((ux-lx)>0 && (uy-ly)>0 && "Invalid indices for mask matrix");
      }
    MaskMatrix ( int lx, int ux, int ly, int uy, const data_type& dfl ):
      MaskMatrix (lx, ux, ly, uy)
      {
        DynamicMatrix::fill(dfl);
      }
    MaskMatrix ( int l, int u ):
      MaskMatrix(l, u, l, u)
      { }
    MaskMatrix ( int l, int u, const data_type& dfl ):
      MaskMatrix(l, u, l, u, dfl)
      { }
    MaskMatrix ( int lx, int ux, int ly, int uy, const DynamicMatrix& rhs ):
      DynamicMatrix(rhs), low_x(lx), up_x(ux), low_y(ly), up_y(uy)
      {
        assert((ux-lx)>0 && (uy-ly)>0 && "Invalid indices for mask matrix");
        assert(rhs.rows()==(ux-lx) && rhs.cols()==(uy-ly) && "Array dimensions do not match");
      }
    MaskMatrix ( int l, int u, const DynamicMatrix& rhs ):
      MaskMatrix(l, u, l, u, rhs)
      {
        assert(rhs.rows()==(u-l) && rhs.cols()==(u-l) && "Array dimensions do not match");
      }
    ~MaskMatrix() = default;
    // Getters
    inline int get_lx (void) const { return low_x; }
    inline int get_ly (void) const { return low_y; }
    inline int get_ux (void) const { return up_x; }
    inline int get_uy (void) const { return up_y; }
    // Accessing operators
    data_type& operator () (int i, int j)
    {
      return DynamicMatrix::operator()( idx_i(i), idx_j(j) );
    }
    const data_type& operator () (int i, int j) const
    {
      return DynamicMatrix::operator()( idx_i(i), idx_j(j) );
    }
    // Assignment operator in order to fill the matrix
    MaskMatrix& operator = (const data_type& rhs)
    {
      DynamicMatrix::fill(rhs);
      return *this;
    }
    // Assignment operator in order to be able to perform a = b (op) c
    MaskMatrix& operator = (const DynamicMatrix& rhs)
    {
      DynamicMatrix::operator=(rhs);
      return *this;
    }
    // Copy-cast
    template <class rhs_type>
    MaskMatrix<data_type>& copy_cast (const MaskMatrix<rhs_type>& rhs)
    {
      DynamicMatrix::operator=(rhs.template cast<data_type>());
      return *this;
    }
    // Block
    constDynamicBlock submatrix ( int lx, int ux, int ly, int uy ) const
    {
      return DynamicMatrix::block( idx_i(lx), idx_j(ly), ux-lx, uy-ly );
    }
  };

  // ###################################################################### //
  // ##### SLIDE MATRIX ################################################### //
  // ###################################################################### //

  template <class data_type, int StorageOrder = RowMajor>
  class SlideMaskMatrix : public MaskMatrix<data_type, StorageOrder>
  {
  protected:
    typedef MaskMatrix<data_type,StorageOrder> Base;
    typedef typename Base::DynamicMatrix DynamicMatrix;
    int n_cut_x;
    int n_cut_y;
  public:
    SlideMaskMatrix(int nxc, int nyc):
      Base(-nxc,nxc+1,-nyc,nyc+1), n_cut_x(nxc), n_cut_y(nyc)
      { }
    SlideMaskMatrix(int nxc, int nyc, const data_type& dfl):
      Base(-nxc,nxc+1,-nyc,nyc+1,dfl), n_cut_x(nxc), n_cut_y(nyc)
      { }
    // Getters
    inline int get_n_cut_x(void) const { return n_cut_x; }
    inline int get_n_cut_y(void) const { return n_cut_y; }
  };

  // ###################################################################### //
  // ##### HALO MATRIX #################################################### //
  // ###################################################################### //

  static constexpr std::array<int,N_BUF> reflect_map = {4, 5, 6, 7, 0, 1, 2, 3};
  static inline int reflect_idx(const int& r)
  {
    return reflect_map[r];
  }

  template <class data_type, int StorageOrder = RowMajor>
  class HaloMaskMatrix : public MaskMatrix<data_type, StorageOrder>
  {
  protected:
    typedef MaskMatrix<data_type,StorageOrder> Base;
    typedef typename Base::DynamicMatrix DynamicMatrix;
    int n_halo_x;       /*!< Number of cut-off cells in the x_dir */
    int n_halo_y;       /*!< Number of cut-off cells in the y_dir */
    int ux_recv[N_BUF], lx_recv[N_BUF];   /*!< Upper and lower x-indices of receiving blocks */
    int uy_recv[N_BUF], ly_recv[N_BUF];   /*!< Upper and lower y-indices of receiving blocks */
    int ux_send[N_BUF], lx_send[N_BUF];   /*!< Upper and lower x-indices of sending blocks */
    int uy_send[N_BUF], ly_send[N_BUF];   /*!< Upper and lower y-indices of sending blocks */
    void inline set_schema()
    {
      // OUTER BUFFER [0]
      lx_recv[TL] = Base::low_x;              ly_recv[TL] = Base::low_y;
      ux_recv[TL] = Base::low_x+n_halo_x;     uy_recv[TL] = Base::low_y+n_halo_y;
      // OUTER BUFFER [1]
      lx_recv[CL] = Base::low_x+n_halo_x;     ly_recv[CL] = Base::low_y;
      ux_recv[CL] = Base::up_x-n_halo_x;      uy_recv[CL] = Base::low_y+n_halo_y;
      // OUTER BUFFER [2]
      lx_recv[BL] = Base::up_x-n_halo_x;      ly_recv[BL] = Base::low_y;
      ux_recv[BL] = Base::up_x;               uy_recv[BL] = Base::low_y+n_halo_y;
      // OUTER BUFFER [3]
      lx_recv[BC] = Base::up_x-n_halo_x;      ly_recv[BC] = Base::low_y+n_halo_y;
      ux_recv[BC] = Base::up_x;               uy_recv[BC] = Base::up_y-n_halo_y;
      // OUTER BUFFER [4]
      lx_recv[BR] = Base::up_x-n_halo_x;      ly_recv[BR] = Base::up_y-n_halo_y;
      ux_recv[BR] = Base::up_x;               uy_recv[BR] = Base::up_y;
      // OUTER BUFFER [5]
      lx_recv[CR] = Base::low_x+n_halo_x;     ly_recv[CR] = Base::up_y-n_halo_y;
      ux_recv[CR] = Base::up_x-n_halo_x;      uy_recv[CR] = Base::up_y;
      // OUTER BUFFER [6]
      lx_recv[TR] = Base::low_x;              ly_recv[TR] = Base::up_y-n_halo_y;
      ux_recv[TR] = Base::low_x+n_halo_x;     uy_recv[TR] = Base::up_y;
      // OUTER BUFFER [7]
      lx_recv[TC] = Base::low_x;              ly_recv[TC] = Base::low_y+n_halo_y;
      ux_recv[TC] = Base::low_x+n_halo_x;     uy_recv[TC] = Base::up_y-n_halo_y;
      // INNER BUFFER [0]
      lx_send[TL] = Base::low_x+n_halo_x;     ly_send[TL] = Base::low_y+n_halo_y;
      ux_send[TL] = Base::low_x+2*n_halo_x;   uy_send[TL] = Base::low_y+2*n_halo_y;
      // INNER BUFFER [1]
      lx_send[CL] = Base::low_x+n_halo_x;     ly_send[CL] = Base::low_y+n_halo_y;
      ux_send[CL] = Base::up_x-n_halo_x;      uy_send[CL] = Base::low_y+2*n_halo_y;
      // INNER BUFFER [2]
      lx_send[BL] = Base::up_x-2*n_halo_x;    ly_send[BL] = Base::low_y+n_halo_y;
      ux_send[BL] = Base::up_x-n_halo_x;      uy_send[BL] = Base::low_y+2*n_halo_y;
      // INNER BUFFER [3]
      lx_send[BC] = Base::up_x-2*n_halo_x;    ly_send[BC] = Base::low_y+n_halo_y;
      ux_send[BC] = Base::up_x-n_halo_x;      uy_send[BC] = Base::up_y-n_halo_y;
      // INNER BUFFER [4]
      lx_send[BR] = Base::up_x-2*n_halo_x;    ly_send[BR] = Base::up_y-2*n_halo_y;
      ux_send[BR] = Base::up_x-n_halo_x;      uy_send[BR] = Base::up_y-n_halo_y;
      // INNER BUFFER [5]
      lx_send[CR] = Base::low_x+n_halo_x;     ly_send[CR] = Base::up_y-2*n_halo_y;
      ux_send[CR] = Base::up_x-n_halo_x;      uy_send[CR] = Base::up_y-n_halo_y;
      // INNER BUFFER [6]
      lx_send[TR] = Base::low_x+n_halo_x;     ly_send[TR] = Base::up_y-2*n_halo_y;
      ux_send[TR] = Base::low_x+2*n_halo_x;   uy_send[TR] = Base::up_y-n_halo_y;
      // INNER BUFFER [7]
      lx_send[TC] = Base::low_x+n_halo_x;     ly_send[TC] = Base::low_y+n_halo_y;
      ux_send[TC] = Base::low_x+2*n_halo_x;   uy_send[TC] = Base::up_y-n_halo_y;
    }
  public:
    // CONSTRUCTOR
    HaloMaskMatrix( const MaskMatrix<data_type>& data_, int nc ):
      MaskMatrix<data_type>(data_), n_halo_x(nc), n_halo_y(nc)
      {
        set_schema();
      }
    HaloMaskMatrix( const MaskMatrix<data_type>& data_, int nc_x, int nc_y ):
      MaskMatrix<data_type>(data_), n_halo_x(nc_x), n_halo_y(nc_y)
      {
        set_schema();
      }
    ~HaloMaskMatrix() = default;
    // Patch-block functionalities
    template <class rhs_type>
    HaloMaskMatrix<data_type>& copy_patch(const MaskMatrix<rhs_type>& rhs, int px, int py)
    {
      DynamicMatrix::block(Base::idx_i(px), Base::idx_j(py), rhs.rows(), rhs.cols())
        = rhs.template cast<data_type>();
      return *this;
    }
    void set_outer_block ( int r, const DynamicMatrix& rhs )
    {
      DynamicMatrix::block(Base::idx_i(lx_recv[r]), Base::idx_j(ly_recv[r]),
        ux_recv[r]-lx_recv[r], uy_recv[r]-ly_recv[r]) = rhs;
    }
    DynamicMatrix get_inner_block ( int r )
    {
      return DynamicMatrix::block(Base::idx_i(lx_send[r]), Base::idx_j(ly_send[r]),
        ux_send[r]-lx_send[r], uy_send[r]-ly_send[r]);
    }
    // Getters
    int get_n_halo_x(void) const { return n_halo_x; }
    int get_n_halo_y(void) const { return n_halo_y; }
  };

  // ###################################################################### //
  // ##### CONVOLUTIONER ################################################## //
  // ###################################################################### //

  template <class data_type>
  class MatrixConvolutioner
  {
  private:
    MaskMatrix<data_type>& result;
    const SlideMaskMatrix<data_type>& slider;
    const HaloMaskMatrix<data_type>& base;
    // Helpers:
    int n_halo_x;   // Actually used only for consistency check ...
    int n_halo_y;   // Actually used only for consistency check ...
    int n_cut_x;
    int n_cut_y;
    int lx;
    int ly;
    int ux;
    int uy;
    data_type def;            // e.g. const real_number default = 0.0
    data_type temp_result;
  public:
    MatrixConvolutioner(MaskMatrix<data_type>& r, const SlideMaskMatrix<data_type>& s,
      const HaloMaskMatrix<data_type>& b, const data_type& d):
      result(r), slider(s), base(b),
      n_halo_x(b.get_n_halo_x()), n_halo_y(b.get_n_halo_y()),
      n_cut_x(s.get_n_cut_x()), n_cut_y(s.get_n_cut_y()),
      lx(r.get_lx()), ly(r.get_ly()), ux(r.get_ux()), uy(r.get_uy()),
      def(d)
      {
        assert( (n_halo_x >= n_cut_x) && "Cut-off in x direction does not match for the convolutioner" );
        assert( (n_halo_y >= n_cut_y) && "Cut-off in y direction does not match for the convolutioner" );
      }
    void convolute(void)
    {
      for ( int i = lx; i<ux; ++i )
      {
        for ( int j = ly; j<uy; ++j )
          convolute(i,j);
      }
    }
    void convolute(int i, int j)
    {
      // TIME WHICH ONE IS BETTER ...
      /*
      temp_result = def;
      for ( int ii = -n_cut_x; ii<=n_cut_x; ++ii )
      {
        for ( int jj = -n_cut_y; jj<=n_cut_y; ++jj )
        {
          temp_result += slider(ii,jj) * base(i+ii,j+jj);
        }
      }
      result(i,j) = temp_result;
      */
      result(i,j) = def +
        ( slider * base.submatrix(-n_cut_x+i, n_cut_x+i+1, -n_cut_y+j, n_cut_y+j+1) ).sum();
    }
  };

} /* namespace ev_matrix */

#endif /* EV_MASK_MATRIX_HPP */
