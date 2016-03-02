// Concurrent writes are threadsafe, and concurrent reads are thread safe.  But no concurrent reads & writes could cause a race.
#ifndef NSG_threadsafe_multimap_h
#define NSG_threadsafe_multimap_h

#include "nsg_util.h"

#include <atomic>
#include <vector>

#define RESERVE_SIZE 2 // reserve 2 places in memory for each bucket.  

using namespace std;

// Spinlocks used to avoid concurrency issues when writing to the same bucket in the hash table
class SpinLock {
public:
  void lock() {
    while(lck.test_and_set(std::memory_order_acquire)) {}
  }
    
  bool test_and_set() {
    return lck.test_and_set(std::memory_order_acquire);
  }
  
  void unlock() {
    lck.clear(std::memory_order_release);
  }
    
private:
    std::atomic_flag lck = ATOMIC_FLAG_INIT;
};

template<typename Key,typename Value,typename Hash=std::hash<Key>>
class ThreadSafeMultiMap {

 public:  
  typedef std::pair<Key,Value> bucket_value; //items in buckets are key/value pairs
  typedef std::vector<bucket_value> bucket_data; //data in bucket is a list of key/value pairs
  typedef typename vector<bucket_value>::iterator bucket_iterator;
  
  ThreadSafeMultiMap(hash_t num_buckets, Hash const& hasher_=Hash())
      : buckets(num_buckets), hasher(hasher_), n_buckets(num_buckets) {
    for(auto i=0;i<num_buckets;++i) {
      buckets[i].reset(new bucket_type);
      buckets[i]->reserve();
    }
  }
    
  pair<bucket_iterator, bucket_iterator> equal_range(
      Key const& key, Value const& default_value=Value()) const {
    return get_bucket(key).equal_range_multi(key);
  }
  
  void insert(Key const& key,Value const& value) {
    get_bucket(key).insert(key,value);
  }
  
  
  hash_t num_buckets() {
    return n_buckets;
  }
  void clear() {
    clear_buckets(0, n_buckets);
  }
  //clear buckets bucket_index0 through (bucket_index1-1)
  void clear_buckets(std::size_t const bucket_index0, std::size_t const bucket_index1) {
    for(auto i=bucket_index0;i<bucket_index1;++i) {
      if(!(buckets[i]->isempty())){
        buckets[i]->clear_bucket();
        buckets[i]->reserve();
      }
    }
  }

 private:
  class bucket_type {  
   public:
    void insert(Key const& key,Value const& value) {
      spinlock.lock();
      data.push_back(bucket_value(key,value));
      empty_flag = false;
      spinlock.unlock();
    }
    void clear_bucket() {
      data.clear();
      empty_flag = true;
    }
    void reserve() {
      data.reserve(RESERVE_SIZE);
    }
    bool isempty() {
      return empty_flag;
    }
    
    pair<bucket_iterator, bucket_iterator> equal_range_multi(Key const& key) {
      if (!empty_flag) {
        bucket_iterator __i = find_if(
            data.begin(),data.end(),[&](bucket_value const& item){ return item.first==key; });
        bucket_iterator __j = __i;
        if (__i != data.end()) {
          bucket_iterator __e = data.end();
          do {
            ++__j;
          } while (__j != __e && ((*__j).first==key));
        }
        return pair<bucket_iterator, bucket_iterator>(__i, __j);
      } else {
        return pair<bucket_iterator, bucket_iterator>(data.end(),data.end());
      }
    }

   private:
    // vector of bucket_values
    bucket_data data;
    // lock to keep concurrency
    SpinLock spinlock;
    // flag indicating if bucket is empty 
    bool empty_flag;
  };
  
  std::vector<std::unique_ptr<bucket_type> > buckets;
  Hash hasher;
  hash_t n_buckets;
  
  bucket_type& get_bucket(Key const& key) const {
    std::size_t const bucket_index = hasher(key) % n_buckets;//maybe change to num_buckets to increase speed?
    return *buckets[bucket_index];
  }    
};

#endif