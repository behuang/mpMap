#ifndef UNIQUE_HEADER_GUARD
#define UNIQUE_HEADER_GUARD
template <typename T, typename V> struct Unique
{
	V value;
	Unique(const V& value = V())
	:value(value)
	{}
	operator V() const
	{
		return value;
	}
};

struct markerEncoding_imp;
typedef Unique<markerEncoding_imp, int> markerEncoding;

struct markerPatternID_imp;
typedef Unique<markerPatternID_imp, int> markerPatternID;

struct funnelEncoding_imp;
typedef Unique<funnelEncoding_imp, int> funnelEncoding;

struct funnelID_imp;
typedef Unique<funnelID_imp, int> funnelID;
#endif