#ifndef NULL_STREAM_H_
#define NULL_STREAM_H_

class NullBuffer : public std::streambuf {
public:
  int overflow(int c) { return c; }
};

class NullStream : public std::ostream {
public:
  NullStream() : std::ostream(&m_sb) {}
private:
  NullBuffer m_sb;
};

#endif /* NULL_STREAM_H_ */
